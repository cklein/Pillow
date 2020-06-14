#
# The Python Imaging Library.
# $Id$
#
# im.show() drivers
#
# History:
# 2008-04-06 fl   Created
#
# Copyright (c) Secret Labs AB 2008.
#
# See the README file for information on usage and redistribution.
#
import os
import shutil
import subprocess
import sys
import tempfile
from shlex import quote

from PIL import Image

_viewers = []


def register(viewer, order=1):
    """
    The :py:func:`register` function is used to register additional viewers.

    :param viewer: The viewer to be registered.
    :param order:
        A negative integer to prepend this viewer to the list,
        or a positive integer to append it.
    """
    try:
        if issubclass(viewer, Viewer):
            viewer = viewer()
    except TypeError:
        pass  # raised if viewer wasn't a class
    if order > 0:
        _viewers.append(viewer)
    elif order < 0:
        _viewers.insert(0, viewer)


def show(image, title=None, **options):
    r"""
    Display a given image.

    :param image: An image object.
    :param title: Optional title. Not all viewers can display the title.
    :param \**options: Additional viewer options.
    :returns: ``True`` if a suitable viewer was found, ``False`` otherwise.
    """
    for viewer in _viewers:
        if viewer.show(image, title=title, **options):
            return 1
    return 0


class Viewer:
    """Base class for viewers."""

    # main api

    def show(self, image, **options):
        """
        The main function for displaying an image.
        Converts the given image to the target format and displays it.
        """

        # save temporary image to disk
        if not (
            image.mode in ("1", "RGBA") or (self.format == "PNG" and image.mode == "LA")
        ):
            base = Image.getmodebase(image.mode)
            if image.mode != base:
                image = image.convert(base)

        return self.show_image(image, **options)

    # hook methods

    format = None
    """The format to convert the image into."""
    options = {}
    """Additional options used to convert the image."""

    def get_format(self, image):
        """Return format name, or ``None`` to save as PGM/PPM."""
        return self.format

    def get_command(self, file, **options):
        """
        Returns the command used to display the file.
        Not implemented in the base class.
        """
        raise NotImplementedError

    def save_image(self, image):
        """Save to temporary file and return filename."""
        return image._dump(format=self.get_format(image), **self.options)

    def show_image(self, image, **options):
        """Display the given image."""
        return self.show_file(self.save_image(image), **options)

    def show_file(self, file, **options):
        """Display the given file."""
        os.system(self.get_command(file, **options))
        return 1


# --------------------------------------------------------------------


class WindowsViewer(Viewer):
    """The default viewer on Windows is the default system application for PNG files."""

    format = "PNG"
    options = {"compress_level": 1}

    def get_command(self, file, **options):
        return (
            'start "Pillow" /WAIT "%s" '
            "&& ping -n 2 127.0.0.1 >NUL "
            '&& del /f "%s"' % (file, file)
        )


if sys.platform == "win32":
    register(WindowsViewer)


class MacViewer(Viewer):
    """The default viewer on MacOS using ``Preview.app``."""

    format = "PNG"
    options = {"compress_level": 1}

    def get_command(self, file, **options):
        # on darwin open returns immediately resulting in the temp
        # file removal while app is opening
        command = "open -a Preview.app"
        command = "({} {}; sleep 20; rm -f {})&".format(
            command, quote(file), quote(file)
        )
        return command

    def show_file(self, file, **options):
        """Display given file"""
        fd, path = tempfile.mkstemp()
        with os.fdopen(fd, "w") as f:
            f.write(file)
        with open(path, "r") as f:
            subprocess.Popen(
                ["im=$(cat); open -a Preview.app $im; sleep 20; rm -f $im"],
                shell=True,
                stdin=f,
            )
        os.remove(path)
        return 1


if sys.platform == "darwin":
    register(MacViewer)


class UnixViewer(Viewer):
    format = "PNG"
    options = {"compress_level": 1}

    def get_command(self, file, **options):
        command = self.get_command_ex(file, **options)[0]
        return "({} {}; rm -f {})&".format(command, quote(file), quote(file))

    def show_file(self, file, **options):
        """Display given file"""
        fd, path = tempfile.mkstemp()
        with os.fdopen(fd, "w") as f:
            f.write(file)
        with open(path, "r") as f:
            command = self.get_command_ex(file, **options)[0]
            subprocess.Popen(
                ["im=$(cat);" + command + " $im; rm -f $im"], shell=True, stdin=f
            )
        os.remove(path)
        return 1


class DisplayViewer(UnixViewer):
    """The ImageMagick ``display`` command."""

    def get_command_ex(self, file, **options):
        command = executable = "display"
        return command, executable


class EogViewer(UnixViewer):
    """The GNOME Image Viewer ``eog`` command."""

    def get_command_ex(self, file, **options):
        command = executable = "eog"
        return command, executable


class XVViewer(UnixViewer):
    """
    The X Viewer ``xv`` command.
    This viewer supports the ``title`` parameter.
    """

    def get_command_ex(self, file, title=None, **options):
        # note: xv is pretty outdated.  most modern systems have
        # imagemagick's display command instead.
        command = executable = "xv"
        if title:
            command += " -name %s" % quote(title)
        return command, executable


if sys.platform not in ("win32", "darwin"):  # unixoids
    if shutil.which("display"):
        register(DisplayViewer)
    if shutil.which("eog"):
        register(EogViewer)
    if shutil.which("xv"):
        register(XVViewer)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Syntax: python ImageShow.py imagefile [title]")
        sys.exit()

    with Image.open(sys.argv[1]) as im:
        print(show(im, *sys.argv[2:]))
