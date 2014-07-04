from helper import unittest, PillowTestCase, tearDownModule

from PIL import Image
from PIL import ImageDraw
from io import BytesIO
import os

FONT_PATH = "Tests/fonts/FreeMono.ttf"
FONT_SIZE = 20


try:
    from PIL import ImageFont
    ImageFont.core.getfont  # check if freetype is available

    class TestImageFont(PillowTestCase):

        def test_sanity(self):
            self.assertRegexpMatches(
                ImageFont.core.freetype2_version, "\d+\.\d+\.\d+$")

        def test_font_with_name(self):
            ImageFont.truetype(FONT_PATH, FONT_SIZE)
            self._render(FONT_PATH)
            self._clean()

        def _font_as_bytes(self):
            with open(FONT_PATH, 'rb') as f:
                font_bytes = BytesIO(f.read())
            return font_bytes

        def test_font_with_filelike(self):
            ImageFont.truetype(self._font_as_bytes(), FONT_SIZE)
            self._render(self._font_as_bytes())
            # Usage note:  making two fonts from the same buffer fails.
            # shared_bytes = self._font_as_bytes()
            # self._render(shared_bytes)
            # self.assertRaises(Exception, lambda: _render(shared_bytes))
            self._clean()

        def test_font_with_open_file(self):
            with open(FONT_PATH, 'rb') as f:
                self._render(f)
            self._clean()

        def test_font_old_parameters(self):
            self.assert_warning(
                DeprecationWarning,
                lambda: ImageFont.truetype(filename=FONT_PATH, size=FONT_SIZE))

        def _render(self, font):
            txt = "Hello World!"
            ttf = ImageFont.truetype(font, FONT_SIZE)
            w, h = ttf.getsize(txt)
            img = Image.new("RGB", (256, 64), "white")
            d = ImageDraw.Draw(img)
            d.text((10, 10), txt, font=ttf, fill='black')

            img.save('font.png')
            return img

        def _clean(self):
            os.unlink('font.png')

        def test_render_equal(self):
            img_path = self._render(FONT_PATH)
            with open(FONT_PATH, 'rb') as f:
                font_filelike = BytesIO(f.read())
            img_filelike = self._render(font_filelike)

            self.assert_image_equal(img_path, img_filelike)
            self._clean()

        def test_render_multiline(self):
            im = Image.new(mode='RGB', size=(300, 100))
            draw = ImageDraw.Draw(im)
            ttf = ImageFont.truetype(FONT_PATH, FONT_SIZE)
            line_spacing = draw.textsize('A', font=ttf)[1] + 8
            lines = ['hey you', 'you are awesome', 'this looks awkward']
            y = 0
            for line in lines:
                draw.text((0, y), line, font=ttf)
                y += line_spacing

            target = 'Tests/images/multiline_text.png'
            target_img = Image.open(target)

            # some versions of freetype have different horizontal spacing.
            # setting a tight epsilon, I'm showing the original test failure
            # at epsilon = ~38.
            self.assert_image_similar(im, target_img, .5)

        def test_rotated_transposed_font(self):
            img_grey = Image.new("L", (100, 100))
            draw = ImageDraw.Draw(img_grey)
            word = "testing"
            font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

            orientation = Image.ROTATE_90
            transposed_font = ImageFont.TransposedFont(
                font, orientation=orientation)

            # Original font
            draw.setfont(font)
            box_size_a = draw.textsize(word)

            # Rotated font
            draw.setfont(transposed_font)
            box_size_b = draw.textsize(word)

            # Check (w,h) of box a is (h,w) of box b
            self.assertEqual(box_size_a[0], box_size_b[1])
            self.assertEqual(box_size_a[1], box_size_b[0])

        def test_unrotated_transposed_font(self):
            img_grey = Image.new("L", (100, 100))
            draw = ImageDraw.Draw(img_grey)
            word = "testing"
            font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

            orientation = None
            transposed_font = ImageFont.TransposedFont(
                font, orientation=orientation)

            # Original font
            draw.setfont(font)
            box_size_a = draw.textsize(word)

            # Rotated font
            draw.setfont(transposed_font)
            box_size_b = draw.textsize(word)

            # Check boxes a and b are same size
            self.assertEqual(box_size_a, box_size_b)

        def test_free_type_font_get_name(self):
            # Arrange
            font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

            # Act
            name = font.getname()

            # Assert
            self.assertEqual(('FreeMono', 'Regular'), name)

        def test_free_type_font_get_metrics(self):
            # Arrange
            font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

            # Act
            ascent, descent = font.getmetrics()

            # Assert
            self.assertIsInstance(ascent, int)
            self.assertIsInstance(descent, int)
            self.assertEqual((ascent, descent), (16, 4))  # too exact check?

        def test_load_path_not_found(self):
            # Arrange
            filename = "somefilenamethatdoesntexist.ttf"

            # Act/Assert
            self.assertRaises(IOError, lambda: ImageFont.load_path(filename))

        def test_default_font(self):
            # Arrange
            txt = 'This is a "better than nothing" default font.'
            im = Image.new(mode='RGB', size=(300, 100))
            draw = ImageDraw.Draw(im)

            target = 'Tests/images/default_font.png'
            target_img = Image.open(target)

            # Act
            default_font = ImageFont.load_default()
            draw.text((10, 10), txt, font=default_font)

            # Assert
            self.assert_image_equal(im, target_img)


except ImportError:
    class TestImageFont(PillowTestCase):
        def test_skip(self):
            self.skipTest("ImportError")


if __name__ == '__main__':
    unittest.main()

# End of file
