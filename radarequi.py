import pexpect
import re
import argparse
import time
import pdb

from signal import SIGINT
from os import listdir, remove
from os.path import isfile, join, exists
from pexpect import TIMEOUT

R2DEC_SIG = '/* r2dec pseudo code output */'
GHIDRA_SIG = 'pdg @ {}\n'
R2CMD_PROMPT = '[0x[0-9a-fA-F]{8}]>'
#R2CMD_PROMPT = '\x1b[0m \x00\r\x1b[38;2;193;156;0m[0x[0-9a-fA-F]{8}]>\x1b[0m \x1b[?1001r\x1b[?1000l''
EXPECTED = [R2CMD_PROMPT, pexpect.TIMEOUT]

ansi_escape = re.compile(r'''
    \x1b  # ESC
    (?:   # 7-bit C1 Fe (except CSI)
        [@-Z\\-_]
    |     # or [ for CSI, followed by a control sequence
        \[
        [0-?]*  # Parameter bytes
        [ -/]*  # Intermediate bytes
        [@-~]   # Final byte
    )
''', re.VERBOSE)

address = re.compile(r'''
   \x00 
   |
   \[0x[0-9a-f]{8}\]>
''', re.VERBOSE)


def cleanline(line):
    clean = ansi_escape.sub('', line.decode())
    clean = address.sub('', clean)
    # may not need to clean returns?
    clean = clean.replace('\r', '')
    return clean


class RadareQUI:
    def __init__(self, file):
        self._child = pexpect.spawn('r2', ['-e io.cache=true', file], timeout=1)
        self._child.setecho(False)

    def clear_buff(self):
        # if self._child.before:
        #    print("Flushing buffer with {}\n\n and\n\n {}".format(self._child.before, self._child.buffer))
        #    self._child.expect(r'.+', timeout=3)
        self._child.buffer = b''

    def close(self):
        self._child.kill(SIGINT)
        self._child.close()
        time.sleep(1)
        if self._child.isalive():
            self._child.terminate()

    def send_analyze(self):
        self._child.sendline('aaa')
        self._child.expect(EXPECTED)

        # clear buffer of analysis messages
        try:
            while True:
                print(self._child.readline().decode())
        except TIMEOUT:
            self.clear_buff()
            print("Analysis Complete!")

    def get_funs_names(self):
        self._child.sendline('afl')
        self._child.expect(EXPECTED)
        funs = dict()

        # obtain entry without the big garbage before it
        entry = self._child.readline()
        entry = entry.split()
        funs[entry[len(entry)-4].decode()] = entry[len(entry)-1].decode()
        del entry

        try:
            while True:
                func_line = self._child.readline()
                func_line = func_line.split()
                funs[func_line[0].decode()] = func_line[len(func_line)-1].decode()
        except TIMEOUT:
            self.clear_buff()
            return funs

    def get_fun_body_ghidra(self, name):
        self._child.sendline('pdg @ {}'.format(name))
        self._child.expect(EXPECTED)

        body = ''

        try:
            for line in self._child:
                body += cleanline(line)
        except TIMEOUT:
            self.clear_buff()
            # since first 382 chars are also garbage
            return body[body.find(GHIDRA_SIG.format(name))+len(GHIDRA_SIG.format(name)):]

    def get_fun_body_r2dec(self, name):
        self._child.sendline('pdd @ {}'.format(name))
        self._child.expect(EXPECTED)

        body = ''

        try:
            for line in self._child:
                body += cleanline(line)
        except TIMEOUT:
            self.clear_buff()
            return body[body.find(R2DEC_SIG):]  # since first chars are garbage until the signature


def main(args):
    # extract the executables
    executables = [exe for exe in listdir(args.path) if isfile(join(args.path, exe))]
    print('Found these files: {}'.format(executables))

    for file in executables:
        # skip non object files
        if '.o' not in file:
            continue

        radare = RadareQUI(join(args.path, file))
        radare.send_analyze()
        function_dict = radare.get_funs_names()
        fcount = 1
        ftotal = len(function_dict)
        print('Disecting File: {}'.format(file))
        ghidra_out = "Output/{}_ghidra_output.txt".format(file)
        r2dec_out = "Output/{}_r2dec_output.txt".format(file)

        if exists(ghidra_out):
            remove(ghidra_out)

        if exists(r2dec_out):
            remove(r2dec_out)

        for function_name in function_dict.values():
            print('Function {} of {} in {}:{}'.format(fcount, ftotal, file, function_name))
            with open(ghidra_out, "a") as ghidra_code_file:
                function_body = radare.get_fun_body_ghidra(function_name)
                ghidra_code_file.write(function_body)
                ghidra_code_file.write('\n\n')

            with open(r2dec_out, "a") as r2dec_code_file:
                function_body = radare.get_fun_body_r2dec(function_name)
                r2dec_code_file.write(function_body)
                r2dec_code_file.write('\n\n')
            fcount += 1
        radare.close()
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="RadareQUI, the Radare Quick User Interface.")
    parser.add_argument('path', help="File path that contains the executables to process")
    main(parser.parse_args())
