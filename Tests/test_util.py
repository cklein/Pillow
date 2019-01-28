import sys

from helper import unittest, PillowTestCase

from PIL import _util

py36 = sys.version_info.major >= 3 and sys.version_info.minor >= 6


class TestUtil(PillowTestCase):

    def test_is_string_type(self):
        # Arrange
        color = "red"

        # Act
        it_is = _util.isStringType(color)

        # Assert
        self.assertTrue(it_is)

    def test_is_not_string_type(self):
        # Arrange
        color = (255, 0, 0)

        # Act
        it_is_not = _util.isStringType(color)

        # Assert
        self.assertFalse(it_is_not)

    def test_is_path(self):
        # Arrange
        fp = "filename.ext"

        # Act
        it_is = _util.isStringType(fp)

        # Assert
        self.assertTrue(it_is)

    @unittest.skipIf(not py36, 'os.path support for Paths added in 3.6')
    def test_path_obj_is_path(self):
        # Arrange
        from pathlib import Path
        fp = Path('filename.ext')

        # Act
        it_is = _util.isPath(fp)

        # Assert
        self.assertTrue(it_is)

    def test_is_not_path(self):
        # Arrange
        filename = self.tempfile("temp.ext")
        fp = open(filename, 'w').close()

        # Act
        it_is_not = _util.isPath(fp)

        # Assert
        self.assertFalse(it_is_not)

    def test_is_directory(self):
        # Arrange
        directory = "Tests"

        # Act
        it_is = _util.isDirectory(directory)

        # Assert
        self.assertTrue(it_is)

    def test_is_not_directory(self):
        # Arrange
        text = "abc"

        # Act
        it_is_not = _util.isDirectory(text)

        # Assert
        self.assertFalse(it_is_not)

    def test_deferred_error(self):
        # Arrange

        # Act
        thing = _util.deferred_error(ValueError("Some error text"))

        # Assert
        self.assertRaises(ValueError, lambda: thing.some_attr)


if __name__ == '__main__':
    unittest.main()
