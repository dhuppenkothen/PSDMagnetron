from pylab import *
import os

data = loadtxt('../data/test_data3.dat')
posterior_sample = atleast_2d(loadtxt('posterior_sample.txt'))

saveFrames = False # For making movies
if saveFrames:
  os.system('rm Frames/*.png')

ion()
for i in range(0, posterior_sample.shape[0]):
  hold(False)
  plot(data[:,0], data[:,1], linestyle="steps-mid")
  hold(True)
  loglog(data[:,0], posterior_sample[i, -data.shape[0]:], 'r')
  xlabel('Frequency', fontsize=16)
  ylabel('Power', fontsize=16)
  draw()
  if saveFrames:
    savefig('Frames/' + '%0.4d'%(i+1) + '.png', bbox_inches='tight')
    print('Frames/' + '%0.4d'%(i+1) + '.png')

ioff()
show()
