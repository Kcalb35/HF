# 头文件

`#include "you-head-file.h"`

在h文件和cpp文件里，注意不要重复定义数据结构

# 文件读写

有`fstream` `ifstream` `ofstream`

```c++
// #include <fstream>

std::fstream file;
file.open("your-file-path",std::ios::in| std::ios::out);
int i;
file >> i;
file << i<<endl;
file.close();
```

# 关于容器

## vector的遍历

```cpp
std::vector<int> li = {1,2,3}
for (auto &item : li){
	cout<<li<<endl;
}
```

注意 `auto &item`的引用



# `int`与`string`的转换

```cpp
// #include <string>
std::to_string(yourIntVar);
```



# Make

### windows环境配置

将`\MinGW\bin`下的`\MinGW\bin\mingw32-make.exe`重命名为`make.exe`