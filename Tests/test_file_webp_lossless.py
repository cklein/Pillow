import pytest
from PIL import Image

from .helper import PillowTestCase, assert_image_equal, hopper

_webp = pytest.importorskip("PIL._webp", reason="WebP support not installed")


class TestFileWebpLossless(PillowTestCase):
    def setUp(self):
        if _webp.WebPDecoderVersion() < 0x0200:
            self.skipTest("lossless not included")

        self.rgb_mode = "RGB"

    def test_write_lossless_rgb(self):
        temp_file = self.tempfile("temp.webp")

        hopper(self.rgb_mode).save(temp_file, lossless=True)

        with Image.open(temp_file) as image:
            image.load()

            self.assertEqual(image.mode, self.rgb_mode)
            self.assertEqual(image.size, (128, 128))
            self.assertEqual(image.format, "WEBP")
            image.load()
            image.getdata()

            assert_image_equal(image, hopper(self.rgb_mode))
