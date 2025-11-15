yMax=4.6       # Hight of the RMSD-RMSF "Y"-axis in Angstroms
resShift=30    # Shift resid on graph

traj=$1
Tskip=0        # skip picoseconds before analyzis
[[ $(echo "scale=3;0${Tend} > 0"|bc) -eq 0 ]] && Tend=10000000
[[ $(echo "scale=3;0${Tskip} > 0"|bc) -eq 0 ]] && Tskip=0
OutDir=./
indexfile=../../index.ndx
ref_struc=../../MIN/pd1_solv_ions.tpr



grp_list="pocket1 pocket2 pocket3 pocket4 pocket5"



gmx_pth=/home/leontyev/programs/bin/gromacs/gromacs-2024.3/bin


Ngrps=$(echo "$grp_list"| wc -w)
grps=$(echo "$grp_list"|sed 's/[ \t]\+/\\n/g')
 
rm -f ref_struc.tpr
$gmx_pth/gmx_d grompp -zero -maxwarn 10 -o ref_struc.tpr -f ../InP.mdp -n ../../index.ndx -p ../../${lignm}-InP.top  -c $refstr -r $refstr 2>err.txt
if [ ! -r ref_struc.tpr ]; then
  echo "ref_struc.tpr was not created, exit"
fi
for grp in ${grp_list}; do
   echo -e "${grp}\n" | $gmx_pth/gmx_d gyrate -f ${traj} -s ${ref_struc} -b $Tskip  -n ${indexfile} -o rg_${Tskip}ps-end_res_${grp}.xvg
done


outfnm="rg_analysis.out"

UseGrps=${grp_list}
echo -e "Drawing SASA graphs for groups:\n${UseGrps}"


python3 -c "
try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
except:
    print('\033[91m' + 'Cannot import \"matplotlib\" module. TRY TO EXECUTE THE COMMAND: \"\033[92mconda activate\033[91m\"' + '\033[0m')

logFile = open('"$outfnm"', 'a')
print ('\nGyration Analysis for ' + '${traj}', file = logFile)
NumGrps = ${Ngrps}
av, std, max, min,imax,imin,tmax,tmin = [0]*NumGrps, [0]*NumGrps,[0]*NumGrps, [0]*NumGrps,[0]*NumGrps, [0]*NumGrps,[0]*NumGrps, [0]*NumGrps # Statistics

def draw_bar_per_residue(slGroups, tskip, val_title, data_color):
    NumGrps = len(slGroups)
    slGroupsRaw = '$(echo ${grp_list})'.split()
    ilDrawGrpIdxs = range(NumGrps) 
    
    fig_subpl, ax = plt.subplots(nrows=NumGrps, ncols=1, sharex='row',figsize=(8,7))
    fig_subpl.suptitle('Gyration Radius of Pockets', fontsize=20, fontweight='bold')
    fig_subpl.supxlabel('time, (ns)',fontsize=18, fontweight='bold',y=0.002)
    fig_subpl.supylabel('Rg, (\u212B)',fontsize=18, fontweight='bold',x=0.005)
    
    xtick_major_spacing=20
    xtick_minor_spacing=10
    plt.subplots_adjust(left=0.11,
                        bottom=0.12, 
                        right=0.985, 
                        top=0.90, 
                        wspace=1.0, 
                        hspace=0.0)

    fig_file = val_title + '_' + tskip + 'ps-end_res' + '.jpg'
    iGrp = 0
    for drawGrp in slGroups:
        axis = ax[iGrp] if NumGrps > 1 else ax
        xvg_file = val_title + '_' + tskip + 'ps-end_res_' + drawGrp + '.xvg'
        print ('Add plot for Group = ' + drawGrp + ', xvg_file = ' + xvg_file)
        X, Y = [], []
        for line in open(xvg_file, 'r'):
            if (line.startswith(('#','@'))): continue  # skip comments in XVG-file
            values = [float(s) for s in line.split()]
            X.append(values[0]/1000)   # Convert Time units ps->ns
            Y.append(values[1] * 10)
        print('Finish reading XVG ', xvg_file)


        NumData = len(X)
        t_ini = X[0]
        t_fin = X[-1]
        print('NumData = ',NumData)

        import numpy as np
        yarr  = np.asarray(Y)
        rms_value = np.sqrt(np.mean(yarr**2))
        print(f'The total {val_title} of group {drawGrp}: {rms_value:.3f}')



        axis.plot(X, Y, color=data_color,linewidth=1.0, label=slGroupsRaw[ilDrawGrpIdxs[iGrp]])
        axis.set_title(drawGrp, fontsize=14, loc='center', y=0.75, bbox=dict(facecolor='white',alpha=0.85,pad=2,boxstyle='round,pad=0.1',linestyle='none'))
        axis.set_xlim(xmin=t_ini,xmax=t_fin)
        axis.set_ylim(ymin=4.5,ymax=8.0)
        axis.tick_params(axis='both', which='major', labelsize=12)


        yarr = np.asarray(Y)
        av[iGrp]  = np.mean(yarr)
        std[iGrp] = np.std(yarr)
        max[iGrp] = np.max(yarr)
        min[iGrp] = np.min(yarr)
        imax[iGrp] = yarr.argmax()
        imin[iGrp] = yarr.argmin()
        tmax = X[imax[iGrp]]
        tmin = X[imin[iGrp]]
        print ('%s: %4.3f +/- %5.4f, min = %4.3f, max = %4.3f , tmin=%7.3f, tmax=%7.3f' % ( drawGrp, av[iGrp], std[iGrp], min[iGrp], max[iGrp], tmin, tmax ) )
        print ('%s: %4.3f +/- %5.4f, min = %4.3f, max = %4.3f , tmin=%7.3f, tmax=%7.3f' % ( drawGrp, av[iGrp], std[iGrp], min[iGrp], max[iGrp], tmin, tmax ) , file = logFile)

        axis.axhline(y=av[iGrp], color=data_color, linestyle='-', linewidth=0.7, label='mean value')
        if slGroupsRaw[ilDrawGrpIdxs[iGrp]] == 'REF:N20_REF:H30':
            continue

        iGrp = iGrp + 1
    print ('Start Drawing' )

    ax_flat = ax.flat if NumGrps > 1 else [ax]
    for axis in ax_flat:
        axis.label_outer()
    logFile.close()
    print ('===== Done with all Groups =====' )
    
    
    fig_subpl.savefig(fig_file,dpi=200)

slGroups = '$(echo ${UseGrps})'.split()
tskip = '${Tskip}'
draw_bar_per_residue(slGroups, tskip, 'rg', 'blue')

exit()








Data = RawData
print ('slGroups:', slGroups)
print ('ilDrawGrpIdxs:', ilDrawGrpIdxs)

try:
    import numpy as np
except:
    print('\033[91m' + 'Cannot import \"numpy\" module. TRY TO EXECUTE THE COMMAND: \"\033[92mconda activate\033[91m\"' + '\033[0m')


logFile = open('"$outfnm"', 'a')
print ('\nGyration Analysis for ' + '${traj}', file = logFile)
av, std, max, min,imax,imin,tmax,tmin = [0]*NumGrps, [0]*NumGrps,[0]*NumGrps, [0]*NumGrps,[0]*NumGrps, [0]*NumGrps,[0]*NumGrps, [0]*NumGrps # Statistics


"



exit















for (( ilam = 0 ; ilam < $NumTIPoints ; ilam++ ))
do
    angfnm=$OutDir/'g_angle_ov'$ilam'.xvg'
    echo defining cis/trans for $angfnm
    if [ $ilam -ne 0 ]; then
      echo "$(sed '/^[@#].*/d' $angfnm | awk '(NR % '"$NstrideOut"' != 1 && NR==FNR && '"$NstrideOut"' > 1) {next}; (NR==FNR){idx=int(NR/'"$NstrideOut"')+1;if ('"$NstrideOut"' == 1) idx=NR; a[idx]=0; if ($2 > -90 && $2 <= 90) a[idx]=1;next}{OFMT = "%i1"; print $0,a[FNR]}' - $cisfnm)" > $cisfnm
    else
      sed '/^[@#].*/d' $angfnm | awk '(NR % '"$NstrideOut"' == 1 || '"$NstrideOut"' == 1){a=0; if ($2 > -90 && $2 <= 90) a=1;printf "%1i\n", a }' > $cisfnm
    fi
done # END OF LOOP OVER xvg-files
echo Output is saved to $cisfnm
awk 'NR=='"$Nskip"' || '"$Nskip"' == 0 && NR == 1{split($0,b);for (i=1;i<=NF;i++) s[i]=0} NR > '"$Nskip"'{for (i=1;i<=NF;i++) {a[i]+=$i; if (b[i] != $i ) {s[i]+=1;b[i]=$i}}} END{printf "# cis-frequencies ['"$Nskip"' steps skipped]: "; for (i=1;i<=NF;i++) printf a[i] OFS; printf "\n";printf "# ratio[trans/cis] ['"$Nskip"' steps skipped]: "; for (i=1;i<=NF;i++) {if (a[i] > 0) printf "%8.3f",(FNR-'"$Nskip"'-a[i])/a[i]; else printf "%s"," 1/0"}; printf "\n"; printf "# tr/cis-switches ['"$Nskip"' steps skipped]: "; for (i=1;i<=NF;i++) printf s[i] OFS; printf "\n"}' $cisfnm > $sumfnm
awk '{a=0;for (i=1;i<=NF;i++) a+=$i;printf "%i\n", a }' $cisfnm >> $sumfnm




