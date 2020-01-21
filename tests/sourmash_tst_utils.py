"Various utilities used by sourmash tests."

from __future__ import print_function
import contextlib
import sys
import os
import tempfile
import shutil
import subprocess
import collections
import pprint

import pkg_resources
from pkg_resources import Requirement, resource_filename, ResolutionError
import traceback
from io import open  # pylint: disable=redefined-builtin
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


SIG_FILES = [os.path.join('demo', f) for f in (
  "SRR2060939_1.sig", "SRR2060939_2.sig", "SRR2241509_1.sig",
  "SRR2255622_1.sig", "SRR453566_1.sig", "SRR453569_1.sig", "SRR453570_1.sig")
]


def scriptpath(scriptname='sourmash'):
    """Return the path to the scripts, in both dev and install situations."""
    # note - it doesn't matter what the scriptname is here, as long as
    # it's some script present in this version of sourmash.

    path = os.path.join(os.path.dirname(__file__), "../")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    path = os.path.join(os.path.dirname(__file__), "../../EGG-INFO/")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    for path in os.environ['PATH'].split(':'):
        if os.path.exists(os.path.join(path, scriptname)):
            return path


def _runscript(scriptname):
    """Find & run a script with exec (i.e. not via os.system or subprocess)."""
    namespace = {"__name__": "__main__"}
    namespace['sys'] = globals()['sys']

    try:
        pkg_resources.load_entry_point("sourmash", 'console_scripts', scriptname)()
        return 0
    except pkg_resources.ResolutionError:
        pass

    path = scriptpath()

    scriptfile = os.path.join(path, scriptname)
    if os.path.isfile(scriptfile):
        if os.path.isfile(scriptfile):
            exec(  # pylint: disable=exec-used
                compile(open(scriptfile).read(), scriptfile, 'exec'),
                namespace)
            return 0

    return -1


ScriptResults = collections.namedtuple('ScriptResults',
                                       ['status', 'out', 'err'])

def runscript(scriptname, args, **kwargs):
    """Run a Python script using exec().

    Run the given Python script, with the given args, in the given directory,
    using 'exec'.  Mimic proper shell functionality with argv, and capture
    stdout and stderr.

    When using :attr:`fail_ok`=False in tests, specify the expected error.
    """
    __tracebackhide__ = True
    sysargs = [scriptname]
    sysargs.extend(args)

    cwd = os.getcwd()
    in_directory = kwargs.get('in_directory', cwd)
    fail_ok = kwargs.get('fail_ok', False)

    try:
        status = -1
        oldargs = sys.argv
        sys.argv = sysargs

        oldout, olderr = sys.stdout, sys.stderr
        sys.stdout = StringIO()
        sys.stdout.name = "StringIO"
        sys.stderr = StringIO()

        os.chdir(in_directory)

        try:
            print('running:', scriptname, 'in:', in_directory, file=oldout)
            print('arguments', sysargs, file=oldout)

            status = _runscript(scriptname)
        except SystemExit as err:
            status = err.code
        except:  # pylint: disable=bare-except
            traceback.print_exc(file=sys.stderr)
            status = -1
    finally:
        sys.argv = oldargs
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdout, sys.stderr = oldout, olderr

        os.chdir(cwd)

    if status is None:
        status = 0

    if status != 0 and not fail_ok:
        print(out)
        print(err)
        assert False, (status, out, err)

    return ScriptResults(status, out, err)


def get_test_data(filename):
    filepath = None
    try:
        filepath = resource_filename(
            Requirement.parse("sourmash"), "sourmash/sourmash/test-data/"\
                + filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), 'test-data',
                                filename)
    return filepath


class TempDirectory(object):
    def __init__(self):
        self.tempdir = tempfile.mkdtemp(prefix='sourmashtest_')

    def __enter__(self):
        return self.tempdir

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            shutil.rmtree(self.tempdir, ignore_errors=True)
        except OSError:
            pass

        if exc_type:
            return False


class RunnerContext(object):
    """
    I am a RunnerContext object from sourmash_tst_utils.

    I have methods 'run_sourmash' and 'run', which run Python scripts.

    Take a look at my 'location', 'last_command' and 'last_result' attributes!

    You can use the 'output' method to build filenames in my temp directory.
    """
    def __init__(self, location):
        self.location = location
        self.last_command = None
        self.last_result = None

    def run_sourmash(self, *args, **kwargs):
        "Run the sourmash script with the given arguments."
        kwargs['fail_ok'] = True
        if 'in_directory' not in kwargs:
            kwargs['in_directory'] = self.location

        cmdlist = ['sourmash']
        cmdlist.extend(args)
        self.last_command = " ".join(cmdlist)
        self.last_result = runscript('sourmash', args, **kwargs)

        if self.last_result.status:
            raise ValueError(self)

        return self.last_result

    def run(self, scriptname, *args, **kwargs):
        "Run a script with the given arguments."
        if 'in_directory' not in kwargs:
            kwargs['in_directory'] = self.location
        self.last_command = " ".join(args)
        self.last_result = runscript(scriptname, args, **kwargs)
        return self.last_result

    def output(self, path):
        return os.path.join(self.location, path)

    def __str__(self):
        s = ""
        if self.last_command:
            s += "Last command run:\n{}\n".format(repr(self.last_command))
            if self.last_result:
                s += "\nLAST RESULT:\n"
                s += "- exit code: {}\n\n".format(self.last_result.status)
                if self.last_result.out:
                    s += "- stdout:\n---\n{}---\n".format(self.last_result.out)
                else:
                    s += '(no stdout)\n\n'
                if self.last_result.err:
                    s += "- stderr:\n---\n{}---\n".format(self.last_result.err)
                else:
                    s += '(no stderr)\n'

        return s


def in_tempdir(fn):
    def wrapper(*args, **kwargs):
        with TempDirectory() as location:
            ctxt = RunnerContext(location)
            newargs = [ctxt] + list(args)
            return fn(*newargs, **kwargs)

    return wrapper


def in_thisdir(fn):
    def wrapper(*args, **kwargs):
        ctxt = RunnerContext(os.getcwd())
        newargs = [ctxt] + list(args)
        return fn(*newargs, **kwargs)

    return wrapper


def run_shell_cmd(cmd, fail_ok=False, in_directory=None):
    cwd = os.getcwd()
    if in_directory:
        os.chdir(in_directory)

    print('running: ', cmd)
    try:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        (out, err) = proc.communicate()

        out = out.decode('utf-8')
        err = err.decode('utf-8')

        if proc.returncode != 0 and not fail_ok:
            print('out:', out)
            print('err:', err)
            raise AssertionError("exit code is non zero: %d" % proc.returncode)

        return (proc.returncode, out, err)
    finally:
        os.chdir(cwd)


@contextlib.contextmanager
def working_dir(path):
    """Changes working directory and returns to previous on exit."""
    prev_cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)
