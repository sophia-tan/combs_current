import matplotlib.pyplot as plt

separation_num = 5
skip=10

y_axis = []
for separation_num in range(50):
    exp_val = 1959620 - skip * 2 * 8667 - separation_num * 8667
    y_axis.append(exp_val)

plt.plot([i for i in range(50)], y_axis)
plt.xlim(10,50)
plt.show()
