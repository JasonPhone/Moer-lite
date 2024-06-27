import matplotlib.pyplot as plt
import numpy as np
import random

hit = []
mis = []
with open("./log.txt", "r") as f:
    lines = f.readlines()
    for line in lines:
        rd = list(map(float, line.split()))
        if rd[0] == 0:
            mis.append((rd[1:4:], rd[4:7]))
        elif rd[0] == 1:
            hit.append((rd[1:4:], rd[4:7]))

# random.shuffle(dirs)
hit = hit[:1000]
mis = mis[:1000]

print(hit[:5])

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
# s：marker标记的大小
# c: 颜色  可为单个，可为序列
# depthshade: 是否为散点标记着色以呈现深度外观。对 scatter() 的每次调用都将独立执行其深度着色。
# marker：样式
for ray in hit[:20]:
    o = ray[0]
    d = [a + b for a in ray[0] for b in ray[1]]
    ax.plot3D(
        xs=[o[0], d[0]],
        ys=[o[1], d[1]],
        zs=[o[2], d[2]],
        zdir="z",
        c="black",
        marker=".",
    )
for ray in mis[:20]:
    o = ray[0]
    d = [a + b for a in ray[0] for b in ray[1]]
    ax.plot3D(
        xs=[o[0], d[0]],
        ys=[o[1], d[1]],
        zs=[o[2], d[2]],
        zdir="z",
        c="red",
        marker=".",
    )
plt.show()
