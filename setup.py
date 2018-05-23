from setuptools import setup, find_packages
#from Cython.Build import cythonize


setup(name='cLoops',
    version='0.10',
    author='Yaqiang Cao',
    author_email='caoyaqiang@picb.ac.cn',
    url='https://github.com/nefclan/cLoops',
    description='Loops calling for ChIA-PET,HiChIP and Hi-C data. Can be applied to similar datasets.',
    classifiers=[
        'Environment :: Console',
        'Programming Language :: Python :: 2.7',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    packages=find_packages(exclude=['tests','docs']),
    long_description=open('README.md').read(),
    #setup_requires=["joblib","numpy","seaborn","pandas","scipy","HTSeq"],
    entry_points={
        'console_scripts': [
            'cLoops=cLoops.pipe:main',
                ],
        },
    #temply disable deLoops for furthur development
    scripts = ["scripts/callStripes","scripts/pet2fingerprint","scripts/pet2juice","scripts/pet2washU"],

    )
