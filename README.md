## 激光器脉冲和光谱演化

* 参考了github地址：[SeveNOlogy7/SimMoLFil: A mode-locked fiber laser simulator 锁模光纤激光器仿真 (github.com)](https://github.com/SeveNOlogy7/SimMoLFil)的代码

* 简单的结合主动锁模激光器的调制方式，得到的输出脉冲和光谱演化的方案

* 尚待优化

  ---

  ---

  

### 文件说明

* components_simplify
  * include the code defined the optical components and electrical components.
  * the definition of parameters such as optical fiber and ISO, etc.
* methods_simplify
  * include the evolution methods
  * the main way which I use is the GNSLE referenced from SeveNOlogy7
  * electrical ways and optical ways include  some basic functions of the main code that necessary