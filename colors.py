import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colors
import math

print 'Script started'
fig = plt.figure()
ax = fig.add_subplot(111)

ratio = 1.0 / 3.0
count = math.ceil(math.sqrt(len(colors.cnames)))
x_count = count * ratio
y_count = count / ratio
x = 0
y = 0
w = 1 / x_count
h = 1 / y_count

print '\tStarting loop...'
for c in colors.cnames:
    pos = (x / x_count, y / y_count)
    ax.add_patch(patches.Rectangle(pos, w, h, color=c))
    ax.annotate(c, xy=pos)
    if y >= y_count-1:
        x += 1
        y = 0
    else:
        y += 1
print '\tloop ended'
plt.show()
plt.clear('all')
print 'Script finished successfully'
