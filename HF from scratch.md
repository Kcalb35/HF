从零开始写HF

# 基组的形式？如何表示原子轨道？

有STO函数和GTO函数



GTO需要轨道系数alpha和ijk

### 3-21g举例

```
O    S #1s 轨道系数+组合系数
      0.3220370000E+03       0.5923939339E-01
      0.4843080000E+02       0.3514999608E+00
      0.1042060000E+02       0.7076579210E+00
O    SP #成键2s2p 轨道系数+s的组合系数+p的组合系数
      0.7402940000E+01      -0.4044535832E+00       0.2445861070E+00
      0.1576200000E+01       0.1221561761E+01       0.8539553735E+00
O    SP #不成键的2s2p 轨道系数+s的组合系数+p的组合系数
      0.3736840000E+00       0.1000000000E+01       0.1000000000E+01
```

对基组进行了数据处理，以一个元素举例
```
3 Li 3
1 0 3
16.11957475 0.1543289673
2.936200663 0.5353281423
0.794650487 0.4446345422
2 0 3
0.6362897469 0.6362897469
0.1478600533 0.1478600533
0.0480886784 0.0480886784
2 1 3
-0.09996722919 0.155916275
0.3995128261 0.6076837186
0.7001154689 0.3919573931
```

`3 Li 3` 代表原子序数+原子名称+轨道数量

`1 0 3` 分布代表主量子数+角量子数+组合的高斯函数个数

紧接着n行轨道系数+组合系数



# 如何定义轨道的数据结构

一个GTO轨道有3个幂指数，轨道中心坐标，轨道系数及其在原子轨道表示的组合系数
```cpp
struct GTO
{
    int ang[3];
    double cartesian[3];
    double orbital_exponent;
    double coefficient;
};

```
一个原子轨道有nlm系数，名字，其组成部分的GTOs，以及坐标中心。
```cpp
struct Orbital
{
    // quantum number
    int n;
    // angular number
    int l;
    // magentic number
    int m;

    std::string name;

    // an orbital consist of GTO s
    std::vector<GTO> component;

    double cartesian[3];
};
```
一个原子由原子序数、名字、原子轨道和坐标构成。
```cpp
struct Atom
{
    // atom number
    int n;
    // atom name
    std::string name;
    // an atom consists of orbitals
    std::vector<Orbital> Orbitals;
    // cartesian position
    double cartesian[3];

};
```

