import os
import re

from setuptools import setup


try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
except Exception:  # pragma: no cover
    _bdist_wheel = None


if _bdist_wheel is not None:
    class bdist_wheel(_bdist_wheel):
        def finalize_options(self):
            super().finalize_options()

            build_number = os.environ.get("BFEE_DOCKING_BUILD_NUMBER")
            if build_number:
                if not re.fullmatch(r"[0-9][0-9A-Za-z_.]*", build_number):
                    raise ValueError(
                        "BFEE_DOCKING_BUILD_NUMBER must start with a digit and "
                        "contain only letters, digits, underscores, or dots"
                    )
                self.build_number = build_number

            platform_tag = os.environ.get("BFEE_DOCKING_PLATFORM_TAG")
            if platform_tag:
                self.root_is_pure = False
                self.plat_name = platform_tag

        def get_tag(self):
            platform_tag = os.environ.get("BFEE_DOCKING_PLATFORM_TAG")
            if platform_tag:
                return "py3", "none", platform_tag
            return super().get_tag()


    setup(cmdclass={"bdist_wheel": bdist_wheel})
else:
    setup()
