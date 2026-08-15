#!/usr/bin/env python3

import argparse, re, subprocess, glob
import numpy as np

### Paste your path to the ORCA executables
ORCA_PATH = "/opt/orca_6_1_0_linux_x86-64_shared_openmpi418/orca"

def make_input_file(args, state):
    ### Parsing tempfile name provided by ORCA (somename.extinp.tmp)
    xyzfile = args.xyzfile[0:-15] + '.xyz'
    #
    # Generating input for the ground-state calculation
    #
    if state=='gs':
            if args.conical:
                ### Input template for the ground-state gradient calculation
                gs_template=f'''!RKS {args.functional} EnGrad {args.basis_set} {args.aux_basis_set}
%pal ncpus {args.ncpus} end
*xyzfile {args.charge} 1 {xyzfile}            
             
             '''
                with open(f'0/gs.inp', 'w') as f:
                    f.write(f'{gs_template}')
            else:
                ### Sometimes we need just orbitals    
                gs_template=f'''!UKS {args.functional} {args.basis_set} {args.aux_basis_set}
%pal ncpus {args.ncpus} end
*xyzfile {args.charge} 1 {xyzfile}
     
            '''
                with open(f'0/gs.inp', 'w') as f:
                    f.write(f'{gs_template}')
    #
    # Generating input for the deltaSCF calculation
    #
    elif state=='s1':
        delta_guess = glob.glob('0/delta.gbw')
        if (args.gs_orbs_at_each_step) or (not delta_guess):
            ### Input orbitals - gs.gbw
            delta_template=f'''! {args.functional} EnGrad TightSCF {args.basis_set} {args.aux_basis_set} DELTASCF UHF {args.converger}
! MORead
%moinp "0/gs.gbw"
%pal ncpus {args.ncpus} end
%SCF ALPHACONF 0,1 END
*xyzfile {args.charge} 1 {xyzfile}

            '''            
            with open(f'0/delta.inp', 'w') as f:
                f.write(f'{delta_template}')
        else:
            ### Input orbitals - delta.gbw
            delta_template=f'''! {args.functional} EnGrad TightSCF {args.basis_set} {args.aux_basis_set} DELTASCF UHF {args.converger}
%pal ncpus {args.ncpus} end
%SCF ALPHACONF 0,1 
     DeltaSCFFromGS FALSE 
END
*xyzfile {args.charge} 1 {xyzfile}

            '''            
            with open(f'0/delta.inp', 'w') as f:
                f.write(f'{delta_template}')   
    #                        
    # Generating input for the triplet calculation
    #
    elif state=='t1':
        triplet_template=f'''! UKS {args.functional} EnGrad TightSCF {args.basis_set} {args.aux_basis_set}
%pal ncpus {args.ncpus} end
*xyzfile {args.charge} 3 {xyzfile}

        '''
        with open(f'0/triplet.inp', 'w') as f:
            f.write(f'{triplet_template}')

def gen_engrad_for_orca(args, engrad_basename):
    ### Parsing number of atoms, contaminated S1 energy and gradients
    with open(f'0/delta.engrad') as f:
        text = f.read()
    natoms = re.findall(r' \d+', text)[0].strip()
    energy_delta = re.findall(r'   [-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', text)[0].strip()
    grad_delta = re.findall(r'      [-+\s]?\d\.\d{12}', text)
    grad_delta = np.array(grad_delta, dtype=float)
    ### Parsing T1 energy and gradients, computing correction if do_not_purify=False
    if not args.do_not_purify:
        with open("0/triplet.engrad") as f:
            text = f.read()
        energy_triplet = re.findall(r'   [-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', text)[0].strip()
        grad_triplet= re.findall(r'      [-+\s]?\d\.\d{12}', text)
        grad_triplet = np.array(grad_triplet, dtype=float)
        energy_delta = str(float(energy_delta) * 2 - float(energy_triplet))
        grad_delta = grad_delta * 2 - grad_triplet
    grad_delta = ''.join([f'{x: .16f}\n' for x in grad_delta])
    with open(f'{engrad_basename}_EXT.engrad', 'w') as f:
        f.write(natoms+'\n'+energy_delta+'\n'+grad_delta)

def calc_penalty(grad_output_file, first_grad_file, second_grad_file, sigma = 3.5, alpha = 0.02): # This function is inherited from the ExtOptORCA-OpenQP wrapper
    with open(first_grad_file) as f:
        text = f.read()
    natoms = re.findall(r' \d+', text)[0].strip()
    e1 = float(re.findall(r'   [-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', text)[0].strip())  # First state E
    g1 = re.findall(r'      [-+\s]?\d\.\d{12}', text)						     # First state gradient
    g1 = np.array(g1, dtype=float)
    with open(second_grad_file) as f:
        text = f.read()
    e2 = float(re.findall(r'[-+\s]?\d+\.\d{12}', text)[0].strip()) # Second state E
    g2 = re.findall(r'[-+\s]?\d\.\d{16}', text)					     # Second state gradient
    g2 = np.array(g2, dtype=float)

    Eeff = (e1 + e2) / 2 + sigma * ( (e2 - e1)**2 / (e2 - e1 + alpha) )
    F = ( (g1 + g2) / 2) + 2 * sigma * ( (e2 - e1) / (e2 - e1 + alpha) - 0.5 * ( (e2 - e1) / (e2 - e1 + alpha) )**2 ) * (g2 - g1)
    F_string = ''.join([f'{F: .16f}\n' for F in F.flatten()])
    with open(f'{grad_output_file}', 'w') as f:
        f.write(natoms+'\n'+f'{Eeff}'+'\n'+F_string)
    

def main():
    parser = argparse.ArgumentParser(
                        prog='DeltaSCF optimizer',
                        description='Python interface to ORCA to compute spin-corrected gradients',)
    parser.add_argument('xyzfile', help='Input geometry XYZ file',)
    parser.add_argument('--conical', action = 'store_true', help='Triggers S0/S1 MECI optimization with penalty function',)
    parser.add_argument('--alpha', default = 0.02, type = float, help='Alpha parameter for penalty function (only for CI optimization)',)
    parser.add_argument('--sigma', default = 3.5, type = float, help='Sigma parameter for penalty function (only for CI optimization)',)
    parser.add_argument('--functional', '-f', help='Density functional for SCF and DeltaSCF computations',)
    parser.add_argument('--charge', '-c', type=int, default=0, help='System charge',)
    parser.add_argument('--basis_set', '-b', help='Basis set',)
    parser.add_argument('--ncpus', '-cpu', type=int, default=12, help='Number of cpus allocated for SCF and DeltaSCF computations',)
    parser.add_argument('--aux_basis_set', '-augb', default='', help='Auxilary basis set for RI-MP2',)
    parser.add_argument('--gs_orbs_at_each_step', action = 'store_true', default=False, help='Calculate ground-state orbitals at each step of geometry optimization',)
    parser.add_argument('--do_not_purify', action = 'store_true', default=False, help='Turn off the spin purification procedure and use spin-contaminated gradients instead',)
    parser.add_argument('--converger', type = str, default='', help='Converging algorothm for DeltaSCF',)
    args = parser.parse_args()

    if not args.conical:         
        basename = args.xyzfile[0:-15]
        # Checking if there is a delta.gbw file from the previous step of optimization
        delta_guess = glob.glob('0/delta.gbw') 
        # Doing the gs calculation if delta.gbw is absent or gs_orbs_at_each_step=True
        if (args.gs_orbs_at_each_step) or (not delta_guess):
            make_input_file(args, state='gs')
            subprocess.run(['bash', '-c', f'{ORCA_PATH} 0/gs.inp >& 0/gs.out'], check=True)
        # Then doing deltaSCF and triplet
        make_input_file(args, state='s1')      
        subprocess.run(['bash', '-c', f'{ORCA_PATH} 0/delta.inp >& 0/delta.out'], check=True)
        make_input_file(args, state='t1')
        subprocess.run(['bash', '-c', f'{ORCA_PATH} 0/triplet.inp >& 0/triplet.out'], check=True)
        gen_engrad_for_orca(args, basename)

    else:
        basename = args.xyzfile[0:-15]
        make_input_file(args, state='gs')
        subprocess.run(['bash', '-c', f'{ORCA_PATH} 0/gs.inp >& 0/gs.out'], check=True)
        make_input_file(args, state='s1')
        subprocess.run(['bash', '-c', f'{ORCA_PATH} 0/delta.inp >& 0/delta.out'], check=True)
        # One may want to optimize MECI using spin-contaminated gradients
        if not args.do_not_purify:
            make_input_file(args, state='t1')
            subprocess.run(['bash', '-c', f'{ORCA_PATH} 0/triplet.inp >& 0/triplet.out'], check=True)
        gen_engrad_for_orca(args, basename)
        calc_penalty(f'{basename}_EXT.engrad', '0/gs.engrad', f'{basename}_EXT.engrad', sigma = args.sigma, alpha = args.alpha)
            
if __name__ == '__main__':
    main()

