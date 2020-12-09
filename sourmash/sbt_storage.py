import abc
from io import BytesIO
import os
import shutil
import sys
import tarfile
from tempfile import NamedTemporaryFile
import zipfile
from abc import ABC
from pathlib import Path

class Storage(ABC):

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
        self.close()

    def close(self):
        pass

    def can_open(self, location):
        return False


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
        newpath = path
        fullpath = os.path.join(self.location, self.subdir, path)

        if os.path.exists(fullpath):
            # check for content, if same return path,
            with open(fullpath, 'rb') as f:
                old_content = f.read()
                if old_content == content:
                    return path

            # different content, need to find new path to save
            newpath = None
            n = 0
            while newpath is None:
                testpath = "{}_{}".format(fullpath, n)
                if os.path.exists(testpath):
                    n += 1
                else:
                    # testpath is available, use it as newpath
                    newpath = "{}_{}".format(path, n)

        fullpath = os.path.join(self.location, self.subdir, newpath)
        dirpath = os.path.dirname(fullpath)
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

        with open(fullpath, 'wb') as f:
            f.write(content)

        return newpath

    def load(self, path):
        path = Path(self.location) / self.subdir / path
        return path.read_bytes()


class ZipStorage(Storage):

    def __init__(self, path):
        self.path = os.path.abspath(path)

        dirname = os.path.dirname(self.path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        self.bufferzip = None

        # Turns out we can't delete/modify an entry in a zipfile easily,
        # so we need to check some things:
        if not os.path.exists(self.path):
            # If the file doesn't exist open it in write mode.
            self.zipfile = zipfile.ZipFile(path, mode='w',
                                           compression=zipfile.ZIP_STORED)
        else:
            # If it exists, open it in read mode and prepare a buffer for
            # new/duplicated items. During close() there are checks to see
            # how the original file needs to be updated (append new items,
            # deal with duplicates, and so on)
            self.zipfile = zipfile.ZipFile(path, 'r')
            self.bufferzip = zipfile.ZipFile(BytesIO(), mode="w")

        self.subdir = None
        subdirs = [f for f in self.zipfile.namelist() if f.endswith("/")]
        if len(subdirs) == 1:
            self.subdir = subdirs[0]

    def _save_to_zf(self, zf, path, content):
        # we repeat these steps for self.zipfile and self.bufferzip,
        # so better to have an auxiliary method
        try:
            info = zf.getinfo(path)
        except KeyError:
            # entry not there yet, write a new one
            newpath = path
        else:
            entry_content = zf.read(info)

            if entry_content == content:
                # skip writing
                return path

            # Trying to write new content:
            # create newpath based on path
            newpath = None
            n = 0
            while newpath is None:
                testpath = "{}_{}".format(path, n)
                try:
                    zf.getinfo(testpath)
                except KeyError:
                    # testpath is available, use it as newpath
                    newpath = testpath
                else:
                    n += 1

        zf.writestr(newpath, content)
        return newpath

    def save(self, path, content):
        # First try to save to self.zipfile, if it is not writable
        # or would introduce duplicates then try to save it in the buffer
        try:
            newpath = self._save_to_zf(self.zipfile, path, content)
        except (ValueError, RuntimeError):
            # Can't write in the zipfile, write in buffer instead
            if self.bufferzip:
                newpath = self._save_to_zf(self.bufferzip, path, content)
            else:
                # Throw error, can't write the data
                raise ValueError("can't write data")

        return newpath

    def _load_from_zf(self, zf, path):
        # we repeat these steps for self.zipfile and self.bufferzip,
        # so better to have an auxiliary method
        try:
            return zf.read(path)
        except KeyError:
            path = os.path.join(self.subdir, path)
            return zf.read(path)

    def load(self, path):
        try:
            return self._load_from_zf(self.zipfile, path)
        except KeyError:
            return self._load_from_zf(self.bufferzip, path)

    def init_args(self):
        return {'path': self.path}

    def close(self):
        # This is a bit complicated, but we have to deal with new data
        # (if the original zipfile is read-only) and possible duplicates.

        if self.bufferzip is None:
            # The easy case: just close the zipfile, nothing else to do
            self.zipfile.close()
        else:
            # The complicated one. Need to consider:
            # - Is there data in the buffer?
            # - If there is, is any of it
            #    * duplicated?
            #    * new data?
            buffer_names = set(self.bufferzip.namelist())
            zf_names = set(self.zipfile.namelist())
            if buffer_names:
                new_data = buffer_names - zf_names
                duplicated = buffer_names.intersection(zf_names)

                if duplicated:
                    # bad news, need to create new file...
                    # create a temporary file to write the final version,
                    # which will be copied to the right place later.
                    tempfile = NamedTemporaryFile()
                    final_file = zipfile.ZipFile(tempfile, mode="w")
                    all_data = buffer_names.union(zf_names)

                    for item in all_data:
                        if item in duplicated or item in buffer_names:
                            # we prioritize writing data from the buffer to the
                            # final file
                            final_file.writestr(item, self.bufferzip.read(item))
                        else:
                            # it is only in the zipfile, so write from it
                            final_file.writestr(item, self.zipfile.read(item))

                    # close the files, remove the old one and copy the final
                    # file to the right place.
                    self.zipfile.close()
                    final_file.close()
                    os.unlink(self.path)
                    shutil.move(tempfile.name, self.path)

                elif new_data:
                    # Since there is no duplicated data, we can
                    # reopen self.zipfile in append mode and write the new data
                    self.zipfile.close()
                    zf = zipfile.ZipFile(self.path, mode='a')
                    for item in new_data:
                        zf.writestr(item, self.bufferzip.read(item))
            # finally, close the buffer and release memory
            self.bufferzip.close()

    @staticmethod
    def can_open(location):
        return zipfile.is_zipfile(location)

    def list_sbts(self):
        return [f for f in self.zipfile.namelist() if f.endswith(".sbt.json")]


class IPFSStorage(Storage):

    def __init__(self, pin_on_add=True, **kwargs):
        import ipfshttpclient
        self.ipfs_args = kwargs
        self.pin_on_add = pin_on_add
        self.api = ipfshttpclient.connect(**self.ipfs_args)

    def save(self, path, content):
        new_obj = self.api.add_bytes(content)
        if self.pin_on_add:
            self.api.pin.add(new_obj)
        return new_obj

        # TODO: the above solution is quick and dirty.
        # we actually want something more organized,
        # like putting all the generated objects inside the same dir.
        # Check this call using the files API for an example.
        # api.files_write("/test/file", io.BytesIO(b"hi"), create=True)
        #
        # This is also required to bring the IPFSStorage closer to what the
        # ZipStorage is doing now.

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
        if not isinstance(content, bytes):
            content = bytes(content)
        self.conn.set(path, content)
        return path

    def load(self, path):
        return self.conn.get(path)

    def init_args(self):
        # TODO: do we want to remove stuff like password from here?
        return self.redis_args

    def __exit__(self, type, value, traceback):
        pass
