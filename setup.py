from setuptools import setup

setup(name='PyB3',
        version='0.10',
        description='Python B3 inputters and outputters',
        url='https://bitbucket.di2e.net/users/kerry.wood/repos/pyb3/browse',
        author='Kerry N. Wood',
        author_email='kerry.n.wood@gmail.com',
        keywords="B3 observaton serializer",
        license='MIT',
        packages=['PyB3'],
        # you'll also need PyRSO (but that is also in GitHub)
        install_requires=['astropy','datetime'],
        include_package_data=True,
        zip_safe=False)
