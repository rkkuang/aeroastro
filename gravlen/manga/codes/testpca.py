#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2016-04-10 11:43:54
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$
 
import numpy as np 
import matplotlib.pyplot as plt
 
def PCA(dataMat,topNfeat=5):
#topNfeat=5 默认选择前五个最大的特征值
#减去均值 
    meanVals = np.mean(dataMat,axis = 0)
    dataMean = dataMat - meanVals
#求协方差方阵 
    conMat = dataMean.T.dot(dataMean)
#求特征值和特征向量
    eigVals,eigVects = np.linalg.eig(conMat)  
#对特征值进行排序  
    eigValInd = np.argsort(eigVals)
    #得到的eigValInd是从小到大的排列，对应的原数据中该元素的索引
    #x = np.array([3, 1, 2])
    #np.argsort(x)
    #array([1, 2, 0])
    #从小到大依次是1,2,3,1对应的索引是1,2对应的索引是2,3对应的索引是0
    eigValInd = eigValInd[:-(topNfeat+1):-1]
    #逆序，从最大到最小的前topNfeat个
#除去不需要的特征向量
    redeigVects=eigVects[:,eigValInd]  
#求新的数据矩阵
    lowdataMat = dataMean.dot(redeigVects)
#求从低维还原回来的数据
    condata = (lowdataMat.dot(redeigVects.T)) + meanVals
#输出降完维德数据加均值
    #因为降维后的数据是一维的了，所以只能加上dataMat整体的平均数进行恢复了
    reducedata=lowdataMat+np.mean(dataMat)
    return reducedata,condata
#100个样本
N=100
x=np.linspace(2,4,N)
y=x+2-4
 
x1=x+(np.random.rand(N)-0.5)*1.5
y1=y+(np.random.rand(N)-0.5)*1.5
 
data = np.array([x1,y1])
print(data.T.shape)
a,b=PCA(data.T,1)

plt.plot(x,y,color='g',linestyle='-',marker='',label='ideal')
plt.plot(x1,y1,color='b',linestyle='',marker='.',label='datawithnoise')
plt.plot(b[:,0],b[:,1],color='r',linestyle='',marker='.',label=u'recon curve')
plt.plot(a[:,0],np.zeros(N),color='k',linestyle='',marker='*',label=u'low dim data')
#在x轴上画数据，y轴全设为0，np.zeros(N)
 
plt.legend(loc='upper left')
plt.axis('equal')
 
plt.show()