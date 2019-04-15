#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 17:57:51 2019

@author: Haohe Liu
"""

import csv
import numpy as np

# name of data file
# 数据集名称
def csv_read(fname,line = 0,trans_int = False):
    birth_data = []
    with open(fname) as csvfile:
        csv_reader = csv.reader(csvfile)  # 使用csv.reader读取csvfile中的文件
        birth_header = next(csv_reader)  # 读取第一行每一列的标题
        if(line == 0):
            for row in csv_reader:  # 将csv 文件中的数据保存到birth_data中
                birth_data.append(row)
        else:
            for _ in range(line):
                birth_data.append(next(csv_reader))
                
    if(trans_int == True):
        birth_data = [[int(x) for x in row] for row in birth_data]  # 将数据从string形式转换为float形式
    birth_data = np.array(birth_data)  # 将list数组转化成array数组便于查看数据结构
    birth_header = np.array(birth_header)
#    print(birth_data.shape)  # 利用.shape查看结构。
#    print(birth_header.shape)
    return birth_data,birth_header.shape, birth_data.shape
    
    