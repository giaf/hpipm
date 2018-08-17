#! /bin/bash

#find . -name "*.S" -exec bash -c 'for file; do cat utils/new_lic_c >temp && tail -n +29 $file >> temp && cat temp > $file; done;' {} +
#find . -name "*.c" -exec bash -c 'for file; do cat utils/new_lic_c >temp && tail -n +29 $file >> temp && cat temp > $file; done;' {} +
#find . -name "*.h" -exec bash -c 'for file; do cat utils/new_lic_c >temp && tail -n +29 $file >> temp && cat temp > $file; done;' {} +
#find . -name "Makefile.*" -exec bash -c 'for file; do cat utils/new_lic.makefile >temp && tail -n +29 $file >> temp && cat temp > $file; done;' {} +
