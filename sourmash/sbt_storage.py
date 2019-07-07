from __future__ import print_function, unicode_literals, division

import abc
from io import BytesIO
import os
import tarfile
import zipfile


class Storage(abc.ABCMeta(str('ABC'), (object,), {'__slots__': ()})):
    # this weird baseclass is compatible with Python 2 *and* 3,
    # we can remove once we support only py3.4+

    @abc.abstractmethod
    def save(self, path, content):
        pass

    @abc.abstractmethod
    def load(self, path):
        pass

    def init_args(self):
        return {}

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass


class FSStorage(Storage):

    def __init__(self, location, subdir):
        self.location = location
        self.subdir = subdir

        fullpath = os.path.join(location, subdir)
        if not os.path.exists(fullpath):
            os.makedirs(fullpath)

    def init_args(self):
        return {'path': self.subdir}

    def save(self, path, content):
        "Save a node/leaf."
        with open(os.path.join(self.location, self.subdir, path), 'wb') as f:
            f.write(content)

        return path

    def load(self, path):
        out = BytesIO()
        with open(os.path.join(self.location, self.subdir, path), 'rb') as f:
            out.write(f.read())

        return out.getvalue()


class ZipStorage(Storage):

    def __init__(self, path=None):
        # TODO: leave it open, or close/open every time?

        if path is None:
            # TODO: Open a temporary file?
            pass

        self.path = os.path.abspath(path)

        dirname = os.path.dirname(self.path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        if os.path.exists(self.path):
            self.zipfile = zipfile.ZipFile(path, 'r')
        else:
            self.zipfile = zipfile.ZipFile(path, mode='w',
                                           compression=zipfile.ZIP_DEFLATED)

    def save(self, path, content):
        self.zipfile.writestr(path, content)
        return path

    def load(self, path):
        return self.zipfile.read(path)

    def init_args(self):
        return {'path': self.path}

    def __exit__(self, type, value, traceback):
        self.zipfile.close()


class TarStorage(Storage):

    def __init__(self, path=None):
        # TODO: leave it open, or close/open every time?

        if path is None:
            # TODO: Open a temporary file?
            pass

        self.path = os.path.abspath(path)

        dirname = os.path.dirname(self.path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        if os.path.exists(self.path):
            self.tarfile = tarfile.open(path, 'r')
        else:
            self.tarfile = tarfile.open(path, 'w:gz')

    def save(self, path, content):
        info = tarfile.TarInfo(path)
        info.size = len(content)

        # TODO: check tarfile mode, if read-only reopen as writable
        self.tarfile.addfile(info, BytesIO(content))

        return path

    def load(self, path):
        content = self.tarfile.getmember(path)
        f = self.tarfile.extractfile(content)
        return f.read()

    def init_args(self):
        return {'path': self.path}

    def __exit__(self, type, value, traceback):
        self.tarfile.close()


class IPFSStorage(Storage):

    def __init__(self, pin_on_add=True, **kwargs):
        import ipfshttpclient
        self.ipfs_args = kwargs
        self.pin_on_add = pin_on_add
        self.read_only = False

        if 'preload' in self.ipfs_args:
            del self.ipfs_args['preload']

        try:
            self.api = ipfshttpclient.connect(**self.ipfs_args)
        except ipfsapi.exceptions.ConnectionError:
            self.api = ipfsapi.connect('ipfs.io', 80)
            self.read_only = True

    def save(self, path, content):
        # api.add_bytes(b"Mary had a little lamb")
        if self.read_only:
            raise NotImplementedError('This is a read-only client. '
                                      'Start an IPFS node to be able to save '
                                      'data.')
        new_obj = self.api.add_bytes(content)
        if self.pin_on_add:
            self.api.pin.add(new_obj)
        return new_obj

        # TODO: the above solution is quick and dirty.
        # we actually want something more organized,
        # like putting all the generated objects inside the same dir.
        # Check this call using the files API for an example.
        # api.files_write("/test/file", io.BytesIO(b"hi"), create=True)

    def load(self, path):
        return self.api.cat(path)

    def init_args(self):
        return self.ipfs_args

    def __exit__(self, type, value, traceback):
        # TODO: do nothing for now,
        # but we actually want something more organized,
        # like putting all the generated objects inside the same dir.
        # Use the files API,
        # add files without flush(),
        # and then flush it here?
        pass


class RedisStorage(Storage):

    def __init__(self, **kwargs):
        import redis
        self.redis_args = kwargs
        self.conn = redis.Redis(**self.redis_args)

    def save(self, path, content):
        self.conn.set(path, content)
        return path

    def load(self, path):
        return self.conn.get(path)

    def init_args(self):
        # TODO: do we want to remove stuff like password from here?
        return self.redis_args

    def __exit__(self, type, value, traceback):
        pass
