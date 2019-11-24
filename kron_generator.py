# coding=utf-8
import sys
import os
import numpy as np
import time

def main(argv):

    ScaleFactor = int(argv[1])
    EdgeFactor = int(argv[2])
    # print("输入参数为S=%d,E=%d" % (ScaleFactor, EdgeFactor))

    Dim = pow(2,ScaleFactor)
    Num = EdgeFactor * Dim
    a = 0.57
    b = 0.19
    c = 0.19
    print("将生成顶点数为%d, 边数为%d 的数据集" % (Dim, Num))

    Idx_vec = np.ones((2,Num))

    ab = a+b
    c_norm = c/(1-ab)
    a_norm = a/ab

    for i in range(0,ScaleFactor):
        i_bit = np.random.rand(Num,) > ab
        k_bit = np.random.rand(Num,) > (c_norm*i_bit + a_norm*np.logical_not(i_bit))
        Idx_vec = Idx_vec + pow(2,i) * np.vstack((i_bit, k_bit))
        
    p = np.random.permutation(Dim)
    Idx_vec = p[Idx_vec.astype(int)]

    p = np.random.permutation(Num)
    Idx_vec = Idx_vec[:,p].T

    print("数据生成完毕！正在写入文件")
    t_stamp = int(time.time()*100)%10000
    bin_name = 's%d.e%d.kron.selfgen.bin' % (ScaleFactor, EdgeFactor)
    if os.path.exists(bin_name):
        bin_name = 's%d.e%d.kron.selfgen.bin.%04d' % (ScaleFactor, EdgeFactor, t_stamp)
    f = open(bin_name,'w')
    Idx_vec.astype('uint32').tofile(f)
    f.close()


if __name__ == '__main__':
    main(sys.argv)
