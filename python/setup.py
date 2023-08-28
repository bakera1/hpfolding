from setuptools import setup, find_packages

readme = open('README.md').read()

setup(name='hp',
      version='1.0.0',
      author='Alexander Baker',
      license='MIT',
      setup_requires=['flake8'],
      author_email='baker.alexander@gmail.com',
      description='Package that implements HP MILP for protein folding',
      url='https://github.com/bakera1/CreditDefaultSwapPricer/',
      # ext_package='isda',
      # py_modules=['isda'],
      classifiers=[
          "Development Status :: 4 - Beta",
          "Programming Language :: Python :: 3.11",
          "License :: OSI Approved :: MIT License",
          "Operating System :: POSIX :: Linux",
      ],
      packages=find_packages(include=['hp', 'hp.*']),
      keywords="h milt solve protein biology",
      include_package_data=True)
