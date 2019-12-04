# https://www.jianshu.com/p/9ce282dc3966?tdsourcetag=s_pctim_aiomsg
# # g++ -o  指定生成可执行文件的名称
# g++ two_pointM_lens_inverse_ray_shoot.cpp  `pkg-config opencv --libs --cflags` -o DisplayImage
# # g++ -c 只编译不链接 

main:two_pointM_lens_inverse_ray_shoot.cpp
	g++ two_pointM_lens_inverse_ray_shoot.cpp  `pkg-config opencv4 --libs --cflags` -o main
mainn:npM_lens_inverse_ray_shoot.cpp 
	g++ npM_lens_inverse_ray_shoot.cpp `pkg-config opencv4 --libs --cflags` -o mainn

# test:main.o T1.o
#     g++ -o test main.o T1.o
# #=============================================
# #make会进行自动推导生成 main.o T1.o 可以不需要写这一部分
# main.o:main.cpp T1.h
#     g++ -c main.cpp
# T1.o:T1.cpp T1.h
#     g++ -c T1.cpp
# #=============================================

# #当然也可以一步到位 \ 是换行连接符 便于Makefile的易读,不用都挤在一行
# test2:
#     g++ -o test2 main.cpp \
# T1.cpp
# #=============================================

# clean:
#     rm test main.o T1.o
# say:
#     echo "112232"

# # 读取到第一个targer test，则test是终极目标
# # 生成目标test ，test需要main.o T1.0
# # 查找到存在目标main.o T1.o 先生成 最后生成test

# #clean和say是标签，并不生成“clean”这个文件，这样的target称之为 “伪目标”
# #伪目标的名字不能和文件名重复
# #如： 当前目录下有一个文件/文件夹 名字为clean,运行make clean则会：
# #       make: `clean' is up to date.
# #为了避免这种情况,可以使用一个特殊的标记“.PHONY”来显示地指明一个目标是“伪目标”
# .PHONY: clean
# clean:
#     rm test main.o T1.o
