# StrucGAN
I am a Python package for generating structure of cluster or molecular using Wasserstein Generative Adversarial Network (wGAN), 
in which I use wasserstein distence to measure the similarity of original distribution and the distribution generated by generative network.  
In my framework, the most important step is coding and decoding structures. In the process of structure coding, the cluster (or molecular) is placed in a cubic box with atoms localized at voxels, and 
each atoms are extended to fill other nearby voxels with a gaussian function. 
If atoms localized on the boundary of two voxels, it will be moved to one voxel by a small displacement.



## Usage

### Supposing the environment variable StrucGANPath=${HOME}/StrucGAN

### Step1: Installing libganstruc
Entering Directory ${StrucGANPath}/GANStrucDeal, and running command 
```bash
python setup.py install
```
to install the package in your machine  

### Step2: Structure encoding, and generating struc.h5 file for training wGAN by script ${StrucGANPath}/tools/gan_deal_struc.py
```bash
python gan_deal_struc.py --xyz2hd5 --xyzpath=Your Directory (the path contains structure files with xyz format)
```
and then copy struc.h5 to Directory ${StrucGANPath}/wGAN

### Step3: Training WGAN
in Directory ${StrucGANPath}/wGAN
```bash
nohup python main.py --cuda 
```

### Step4: Generating Structure by wGAN
copying the latest GAN model's parameter file netG_epoch_*.pth in Directory ${StrucGANPath}/samples to Directory ${StrucGANPath}/sampling 
and entering ${StrucGANPath}/sampling  
```bash
python main.py --netG=netG_epoch_*.pth
```

### Step5: Decoding structures form struc.h5 file
in Directory ${StrucGANPath}/sampling
```bash
python gan_deal_struc.py --hd52xyz
```
