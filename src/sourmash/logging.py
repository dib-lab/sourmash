import sys
from io import StringIO

_quiet = False
_debug = False
def set_quiet(val, print_debug=False):
    global _quiet, _debug
    _quiet = bool(val)
    _debug = bool(print_debug)


def print_results(s, *args, **kwargs):
    if _quiet:
        return

    print(s.format(*args, **kwargs), file=sys.stdout)
    sys.stdout.flush()


def notify(s, *args, **kwargs):
    "A simple logging function => stderr."
    if _quiet:
        return

    print(u'\r\033[K', end=u'', file=sys.stderr)
    print(s.format(*args, **kwargs), file=sys.stderr,
          end=kwargs.get('end', u'\n'))
    if kwargs.get('flush'):
        sys.stderr.flush()


def short_notify(msg_template, *args, end='\n', screen_width=80, **kwargs):
    """Shorten the notification message so that it fits on one line.

    Good for repeating notifications with end='\r' especially...
    """
    w = screen_width
    if args or kwargs:
        msg = msg_template.format(*args, **kwargs)
    else:
        msg = msg_template

    if len(msg) > w:
        truncate_len = len(msg) - w + 3
        msg = '<<<' + msg[truncate_len:]

    notify(msg, end=end)


def debug(s, *args, **kwargs):
    "A debug logging function => stderr."
    if _quiet or not _debug:
        return

    print(u'\r\033[K', end=u'', file=sys.stderr)
    print(s.format(*args, **kwargs), file=sys.stderr,
          end=kwargs.get('end', u'\n'))
    if kwargs.get('flush'):
        sys.stderr.flush()


def error(s, *args, **kwargs):
    "A simple error logging function => stderr."
    print(u'\r\033[K', end=u'', file=sys.stderr)
    print(s.format(*args, **kwargs), file=sys.stderr)
    if kwargs.get('flush'):
        sys.stderr.flush()


def test_notify():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        notify(u'hello, world')
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world\n' in saveerr.getvalue()


def test_notify_flush():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        notify(u'hello, world', flush=True)
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world' in saveerr.getvalue()


def test_notify_end():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        notify(u'hello, world', end=u'FOO')
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, worldFOO' in saveerr.getvalue()


def test_notify_quiet():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = True
        notify(u'hello, world')
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world' not in saveerr.getvalue()


def test_error():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        error(u'hello, world')
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world\n' in saveerr.getvalue()


def test_error_flush():
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = False
        error(u'hello, world', flush=True)
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world' in saveerr.getvalue()


def test_error_quiet():
    # error should still output even if _quiet is True
    global _quiet

    qsave = _quiet
    saveerr, sys.stderr = sys.stderr, StringIO()
    try:
        _quiet = True
        error(u'hello, world')
    finally:
        _quiet = qsave
        saveerr, sys.stderr = sys.stderr, saveerr

    print(type(saveerr))
    assert 'hello, world' in saveerr.getvalue()
