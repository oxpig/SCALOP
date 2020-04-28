from distutils.core import setup
setup(name='scalop',
      version='1.0',
      description='Sequence-based antibody CAnonical LOoP structure annotation',
      author='Wing Ki Wong',
      author_email='opig@stats.ox.ac.uk',
      url='http://opig.stats.ox.ac.uk/webapps/SCALOP', 
      packages=['scalop',
      			 'scalop.anarci',
      			 'scalop.prosci',
      			 'scalop.prosci.util',
      			 'scalop.prosci.loops'],
      package_dir={'scalop': 'lib/python/scalop',
      				 'scalop.anarci': 'lib/python/scalop/anarci',
      				 'scalop.prosci': 'lib/python/scalop/prosci',
      				 'scalop.prosci.util': 'lib/python/scalop/prosci/util',
      				 'scalop.prosci.loops': 'lib/python/scalop/prosci/loops'},
      package_data={
          'scalop': ['database/*'],
          'scalop.anarci': ['dat/HMMs/ALL*']
      },
      scripts=['bin/SCALOP'],
      license="BSD 3-clause"
     )
