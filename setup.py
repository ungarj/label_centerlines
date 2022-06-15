from setuptools import setup

# get version number
# from https://github.com/mapbox/rasterio/blob/master/setup.py#L55
with open("label_centerlines/__init__.py") as f:
    for line in f:
        if line.find("__version__") >= 0:
            version = line.split("=")[1].strip()
            version = version.strip('"')
            version = version.strip("'")
            continue

setup(
    name="label_centerlines",
    version=version,
    description="extract centerlines from polygons",
    author="Joachim Ungar",
    author_email="joachim.ungar@gmail.com",
    url="https://github.com/ungarj/label_centerlines",
    license="MIT",
    packages=["label_centerlines"],
    entry_points={"console_scripts": ["label_centerlines=label_centerlines.cli:main"]},
    install_requires=[
        "click",
        "Fiona>=1.7.0",
        "networkx>=2.1",
        "scipy>=0.17",
        "Shapely>=1.5",
        "tqdm",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: GIS",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6",
    ],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
)
