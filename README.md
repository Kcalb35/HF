# A simple achieve for Hatree-Fock method in c++
## how to build
```bash
git clone https://github.com/Kcalb35/HF
cd HF

# clone googletest
mkdir lib
cd lib
git clone https://github.com/google/googletest.git
cd ..

# build
mkdir build
cd build
cmake .
make
```

## todo
- [x] parse d type atom orbital
- [ ] optimize output
- [x] optimize calculation
- [ ] optimize command line params