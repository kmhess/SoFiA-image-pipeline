import os
import sys
import unittest

sys.path.append('./src')
from src import image_pipeline


class TestRunPipeline(unittest.TestCase):
    def setUp(self):
        """Verify test catalogue exists

        """
        self.assertTrue(
            os.path.exists('tests/data/UGC7012_cat.xml'),
            "Test catalogue does not exist at tests/data/"
        )

    def test_run_image_pipeline_on_catalogue(self):
        """Run the pipeline on provided test data.

        """
        image_pipeline.main(["-c", "tests/data/UGC7012_cat.xml"])


if __name__ == '__main__':
    unittest.main()
