import numpy as np
import sys
import glob
import os

try:
    from  libganstruc import  encode_struc_3d, decode_struc_3d
except:
    print 'Please install libganstruc'
    sys.exit(0)

try:
    import h5py
except:
    print 'Please install h5py'
    sys.exit(0)


class  GANStruc(object):
    def __init__(self):
        pass


    def read_struc_xyz(self,filename):
        f = open(filename,'r')
        aa = f.readlines()
        natom = int(aa[0])
        element = []
        xyz = []
        ene = 0.0
        for  i in range(natom):
            element.append(aa[ i + 2 ].split()[0])
            xyz.append(aa[ i + 2 ].split()[1:])
        xyz = np.array(xyz,float)
        return natom,element,ene,xyz
    
    def get_volumetric_from_xyz(self, natom, element, xyz, num_grid = 32):
        ele = {'H':1,'He':2,'Li':3,'B':1,'C':6,'N':7,'O':8,'F':9}
        step = 0.3
        bb = encode_struc_3d(num_grid, xyz, step, 0.5)
        return  bb.tolist()
    
    def write_hd5_files(self,path):
        struc = []
        energy = []
        xyzfile = glob.glob(path + '/*.xyz')
        #print xyzfile
        #xyzfile = sorted(xyzfile, key=lambda x: int(x.split('.')[0]))
        for xyz in xyzfile:
            stype = 2
            natom,element,ene,xyz = self.read_struc_xyz(xyz)
            aa = self.get_volumetric_from_xyz(natom,element,xyz)
            energy.append(stype)
            struc.append(aa)
        struc = np.array(struc,float)
        
        f1 = h5py.File('struc.h5','w')	
        f2 = h5py.File('struc_label.h5','w')	
        f1['data'] = struc
        f1['label'] = energy          
        #f2['label'] = energy          
        f1.close()
        f2.close()
    
    
    def read_hd5_files_2D(self,finput):
        if not os.path.exists('./struc_xyz'):
            os.mkdir('struc_xyz')
        f = h5py.File(finput,'r')
        a = f['data'][:]
        f.close()
        (n1, n2, n3, n4) = a.shape
        for  ii  in range(n1):
       	    pos = []
            temp = a[ii,0,:,:]
            self.print_array(ii,temp)
            b=np.where((temp >  0.9))
            for  i in range(len(b[0])):
                pos.append((b[0][i],b[1][i],0.0))
                pos = np.array(pos,int)
                pos = pos*0.3
                foutput = 'struc_xyz/gan_struc_' + str(ii + 1) + '.xyz'
                f = open(foutput,'w')
                f.write(str(len(b[0])) + '\n')
                f.write('\n')
                for  i in range(len(b[0])):
                    f.write('B  ')
                    f.write('%15.11f %15.11f %15.11f\n' % tuple(pos[i]))
                f.close()
    
    def print_array2D(self,n,x):
        f = open('struc_dat/structure_%d.dat' %  (n + 1),'w')
        (n,m) = x.shape
        for i in range(n):
	    for j in range(m):
                f.write(('%5.3f  ' % x[i,j]))
        f.write('\n')
        f.close()
    
    def print_array3D(self, ii, x):
        (n,m,k) = x.shape
        for i in range(k):
    	    f = open('struc_dat/structure_%d_%d.dat' %  (ii, i + 1),'w')
            for i1 in range(n):
                for i2 in range(m):
                    f.write(('%5.3f  ' % x[i1,i2,i]))
                f.write('\n')
            f.close()
    
    def read_hd5_files_3D(self, finput):
        if not os.path.exists('./struc_xyz'):
            os.mkdir('struc_xyz')
        print finput
        f = h5py.File(finput,'r')
        a = f['data'][:]
        f.close()
        discut = 0.3
        na = 80
        print a.shape
        (n1, n2, n3, n4, n5) = a.shape
        for  ii  in range(n1):
            temp = a[ii,0,:,:,:]
            #self.print_array3D(ii,temp)
            (nas, pos) = decode_struc_3d(discut, na, temp)
            if nas == na:
                foutput = 'struc_xyz/gan_struc_' + str(ii + 1) + '.xyz'
                f = open(foutput,'w')
                f.write(str(na) + '\n')
                f.write('\n')
                for  i in range(na):
                    f.write('B  ')
                    f.write('%15.11f %15.11f %15.11f\n' % tuple(pos[i]))
                f.close()
                #sys.exit(0)
	
