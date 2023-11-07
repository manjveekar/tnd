#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: manjveekar
"""


import subprocess


def atom_contacts(coor,resi):
    contacts=0
    for i in coor:
        for j in resi:
            dist=((i[0]-j[0])**2+(i[1]-j[1])**2+(i[2]-j[2])**2)**0.5
            if dist<4.5:
                contacts+=1
    return contacts
def find_edges(pdb):
    aa,res_name,net_file={},{},open('resinet/'+pdb,'w')
    for i in open('PDB/'+pdb+'.pdb','r'):
        if i[22:27] not in aa:
            aa[i[22:27]]=[]
            res_name[i[22:27]]=i[17:20]
        aa[i[22:27]].append([float(i[30:38]),float(i[38:46]),float(i[46:54])])
    nums=sorted(aa.keys())
    for i in range(len(nums)):
        for j in range(len(nums)-i-1):

            if i+1==i+j+1:
                contacts=0
            else:
                contacts=atom_contacts(aa[nums[i]],aa[nums[i+j+1]])
            net_file.write(nums[i]+' '+nums[i+j+1]+' '+res_name[nums[i]]+' '+res_name[nums[i+j+1]]+' '+str(contacts)+'\n')
    net_file.close()
def resinet():
    subprocess.call(['mkdir','resinet'])
    for i in open('input_pdb','r'):
        find_edges(i.split()[0])
#resinet()

def store_edge_weight():
    subprocess.call(['mkdir','network'])
    li,max_con=[],{}
    for i in open('max_contacts','r'):
        max_con[i.split()[0]]=int(i.split()[1])
    for i in open('input_pdb','r'):
        li.append(i.split()[0])
    for i in li:
        out=open('network/'+i,'w')
        for j in open('resinet/'+i,'r'):
            p=max_con[j.split()[2]+'_'+j.split()[3]]
            r1,r2=i[-1]+'_'+j.split()[2]+'_'+j.split()[0].replace(' ',''),i[-1]+'_'+j.split()[3]+'_'+j.split()[1].replace(' ','')
            r3=j.split()[4]
            out.write(r1+' '+r2+' '+r3+' '+str(p)+' '+str(float(j.split()[4])/p)+'\n')
        out.close()
#store_edge_weight()

#get tma
def get_alignments():
    open('tma_output','w')
    pdbs=[]
    for i in open('input_pdb','r'):
        pdbs.append(i[:-1])
    pairs=open('pairs','w')
    for j in range(len(pdbs)):
        for k in range(len(pdbs)-1-j):
            pdb1=pdbs[j]
            pdb2=pdbs[j+k+1]
            pairs.write(pdb1+' '+pdb2+'\n')
            out=open('tma_output','a')
            subprocess.call(['./tma','PDB/'+pdb1+'.pdb','PDB/'+pdb2+'.pdb'],stdout=out)
            out.close()
    pairs.close()
#get_alignments()

def get_seq(pdb1,pdb2):
    seq1,seq2,sequence=[],[],{}
    for i in open('PDB/'+pdb1+'.pdb','r'):
        temp=i[21]+'_'+i[17:20]+'_'+(i[22:27]).replace(' ','')
        if i[12:16]==' CA ' and  temp not in seq1:
            seq1.append(temp)
    for i in open('PDB/'+pdb2+'.pdb','r'):
        temp=i[21]+'_'+i[17:20]+'_'+(i[22:27]).replace(' ','')
        if i[12:16]==' CA ' and  temp not in seq2:
            seq2.append(temp)
    for i in open('alignments/'+pdb1+'_'+pdb2,'r'):
        if i.split()[0]=='-':
            seq2.pop(0)
        elif i.split()[1]=='-':
            seq1.pop(0)
        elif i.split()[2]=='ALIGNED':
            sequence[seq1[0]]=seq2[0]
            seq1.pop(0)
            seq2.pop(0)
    return sequence  
def topo_eq():
    subprocess.call(['mkdir','alignments'])
    rmsd=open('tma_output','r').readlines()
    count=0
    pdbs=[]
    for i in open('input_pdb','r'):
        pdbs.append(i[:-1])
    for j in range(len(pdbs)):
        for k in range(len(pdbs)-1-j):
            out=open('alignments/'+pdbs[j]+'_'+pdbs[j+k+1],'w')
            seq1=(rmsd[(count*26)+22])
            align=(rmsd[(count*26)+23])
            seq2=(rmsd[(count*26)+24])
            for i in range(len(align)):#gives wether aligned or gap
                if align[i]==':' or align[i]=='.':
                    out.write(seq1[i]+' '+seq2[i]+' '+'ALIGNED\n')
                elif seq1[i]=='-':
                    out.write(seq1[i]+' '+seq2[i]+' '+'SEQ\n')
                elif seq2[i]=='-':
                    out.write(seq1[i]+' '+seq2[i]+' '+'SEQ2\n')
            out.close()
            count+=1
            seq=get_seq(pdbs[j],pdbs[j+k+1])
            out=open('alignments/'+pdbs[j]+'_'+pdbs[j+k+1]+'_aligned','w')
            for i in seq:
                out.write(i+' '+seq[i]+' ALI\n')#gives all aligned residues
            out.close()
#topo_eq()
    
def get_norm(pdb,pos):
    count,fgap,fall,gaps,seq=0,0,0,[],[]
    for i in open('PDB/'+pdb+'.pdb','r'):
        res=i[21]+'_'+i[17:20]+'_'+(i[22:27]).replace(' ','')
        if i[12:16]==' CA ' and  res not in seq:
            count+=1
            seq.append(res)
            if count in pos:
                gaps.append(res)
    for i in open('network/'+pdb,'r'):
        if i.split()[1] in gaps or i.split()[0] in gaps:
            fgap+=(float(i.split()[2]))**2
        fall+=(float(i.split()[2]))**2
    nd=fgap**0.5/fall**0.5
    return nd
def compute_tnd():
    line={}
    for i in open('legend','r'):
        temp=i.split()[3]+'_'+i.split()[4]
        line[temp]=i.split()[0]+' '+i.split()[1]+' '+i.split()[2]
    out=open('output_tnd','w')
    tinfo=open('output_tnd_params','w')
    ndsa={}
    for i in open('input_nds','r'):
        ndsa[i.split()[0]+'_'+i.split()[1]]=float(i.split()[3])
    for i in open('pairs','r'):
        pdb1=i.split()[0]
        pdb2=i.split()[1]
        g1,g2,count,pos1,pos2=0,0,0,[],[]
        for j in open('alignments/'+pdb1+'_'+pdb2,'r'):
            count+=1
            if j.split()[2]=='SEQ':
                g1+=1# a gap, corresponding residue in other sequecne is a gap residue
                pos2.append(count-g2)
            elif j.split()[2]=='SEQ2':
                g2+=1
                pos1.append(count-g1)
        nds=ndsa[pdb1+'_'+pdb2]
        nl1=get_norm(pdb1,pos1)
        nl2=get_norm(pdb2,pos2)
        nlcom=float(count-g1-g2)/((count-g1)*(count-g2))**0.5
        tnd=(nds*(1+nl1)*(1+nl2))/(nlcom)
        out.write(line[pdb1+'_'+pdb2]+' '+i[:-1]+' '+str(tnd)+'\n')
        tinfo.write(line[pdb1+'_'+pdb2]+' '+pdb1+' '+pdb2+' '+' '.join([str(x) for x in [tnd,nl1,nl2,nlcom]])+'\n')
        #break
    out.close()
#compute_tnd()
