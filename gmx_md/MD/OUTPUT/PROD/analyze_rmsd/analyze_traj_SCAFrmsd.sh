traj=$1
Tskip=0
[[ $(echo "scale=3;0${Tend} > 0"|bc) -eq 0 ]] && Tend=10000000
[[ $(echo "scale=3;0${Tskip} > 0"|bc) -eq 0 ]] && Tskip=0
OutDir=./
indexfile=../../index.ndx


grp_list="pocket1 pocket2 pocket3 pocket4 pocket5 Bbone"

gmx_pth=/home/leontyev/programs/bin/gromacs/gromacs-2024.3/bin

Ngrps=$(echo "$grp_list"| wc -w)
grps=$(echo "$grp_list"|sed 's/[ \t]\+/\\n/g')
 

echo -e "Bbone\n$grps\n" |$gmx_pth/gmx_d rms -pbc -fit rot+trans -b ${Tskip} -e ${Tend} -ng $Ngrps -s ../../MIN/pd1_solv_ions.tpr -f $traj -n ${indexfile} -o $OutDir/rmsd.xvg
cd $OutDir
sumfnm=SCAF_residues_rmsd.out
Nskip=0    # Skip Nskip Output data points before start averaging number
rmsd_row=$(sed '/^[@#].*/d' rmsd.xvg | awk 'NR=='"$Nskip"' || '"$Nskip"' == 0 && NR == 1{ndata=1;for (i=2;i<=NF;i++) s[i]=$i^2} NR > '"$Nskip"'{ndata+=1;for (i=2;i<=NF;i++) {s[i]+=$i^2}} END{for (i=2;i<=NF;i++) printf "%5.3f ",      sqrt(s[i]/ndata) OFS; printf "\n"}' -)
rmsdA_row=$(sed '/^[@#].*/d' rmsd.xvg | awk 'NR=='"$Nskip"' || '"$Nskip"' == 0 && NR == 1{ndata=1;for (i=2;i<=NF;i++) s[i]=$i^2} NR > '"$Nskip"'{ndata+=1;for (i=2;i<=NF;i++) {s[i]+=$i^2}} END{for (i=2;i<=NF;i++) printf "%5.2f ",10.0*sqrt(s[i]/ndata) OFS; printf "\n"}' -)
echo rmsd_row=$rmsd_row
rmsd=( $rmsd_row )
rmsdA=( $rmsdA_row )
legend=( $(sed -n 's/^[@][ \t]\+s[0-9]\+[ \t]\+legend[ \t]\+\"\([^ \t\"]\+\)\".*/\1/p' rmsd.xvg) )
echo legend=${legend[@]}
echo  -e "\nGRP_NAME\t\t\tRMSD_NM\t\t\tRMSD_A"
echo  -e   "GRP_NAME\t\t\tRMSD_NM\t\t\tRMSD_A" > $sumfnm
for (( ig = 0 ; ig < ${#legend[@]} ; ig++ )); do
  echo -e "${legend[ig]}\t\t\t${rmsd[ig]}\t\t\t${rmsdA[ig]}"
  echo -e "${legend[ig]}\t\t\t${rmsd[ig]}\t\t\t${rmsdA[ig]}" >> $sumfnm
done


  dynfnm="rmsd.xvg"
  dynfig="rmsd.png"
  outfnm="rmsd_analysis.out"
  echo "Plot RMSD Dynamics for ${traj}"
  DrawGrps="pocket1 pocket2 pocket3 pocket4 pocket5 Bbone"

GrpNames=$(sed -n 's/^[ \t]*\[[ \t]*\([^ \t]\+\)[ \t]*\].*/\1/p' ${indexfile})
echo -e "All analyzed groups:\n${GrpNames}"

UseGrps=''
DrawGrpIdxs=''
idx=0
for DrawGrpNm in $DrawGrps; do
   idx=0
   for GrpNm in $grp_list; do
     if [ "${GrpNm}" == "${DrawGrpNm}" ];then
       UseGrps="${UseGrps}\n${GrpNm}"
       DrawGrpIdxs="${DrawGrpIdxs}\n${idx}"
       echo -e "Use group ${GrpNm}"
     fi
     idx=`expr $idx + 1`
done
done
echo -e "Drawing RMSD graphs for groups:\n${UseGrps}"


python3 -c "
try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
except:
    print('\033[91m' + 'Cannot import \"matplotlib\" module. TRY TO EXECUTE THE COMMAND: \"\033[92mconda activate\033[91m\"' + '\033[0m')

slGroupsRaw = '$(echo ${grp_list})'.split()
slGroups = '$(echo ${UseGrps})'.split()
NumGrps = len(slGroups)
ilDrawGrpIdxs = [int(s) for s in '$(echo ${DrawGrpIdxs})'.split()] 

RawData, X =[], []
for line in open('${dynfnm}', 'r'):
    if (line.startswith(('#','@'))): continue  # skip comments in XVG-file
    values = [float(s) for s in line.split()]
    RawData.append(values)
    X.append(values[0]/1000)   # Convert Time units ps->ns
NumData = len(X)
t_ini = X[0]
t_fin = X[NumData-1]
print('NumData = ',NumData)


Data = RawData
print ('slGroups:', slGroups)
print ('ilDrawGrpIdxs:', ilDrawGrpIdxs)

try:
    import numpy as np
except:
    print('\033[91m' + 'Cannot import \"numpy\" module. TRY TO EXECUTE THE COMMAND: \"\033[92mconda activate\033[91m\"' + '\033[0m')


logFile = open('"$outfnm"', 'a')
print ('\nRMSD Analysis for ' + '${traj}', file = logFile)
av, std, max, min,imax,imin,tmax,tmin = [0]*NumGrps, [0]*NumGrps,[0]*NumGrps, [0]*NumGrps,[0]*NumGrps, [0]*NumGrps,[0]*NumGrps, [0]*NumGrps # Statistics

fig_subpl, axis = plt.subplots(nrows=NumGrps, ncols=1, sharex='row',figsize=(6,5))
fig_subpl.suptitle('RMSD Dynamics', fontsize=20, fontweight='bold')
fig_subpl.supxlabel('time, (ns)',fontsize=18, fontweight='bold',y=0.002)
fig_subpl.supylabel('RMSD, (\u212B)',fontsize=18, fontweight='bold',x=0.005)
plt.xlim(t_ini, t_fin)
xtick_major_spacing=20
xtick_minor_spacing=10
plt.subplots_adjust(left=0.11,
                    bottom=0.12, 
                    right=0.985, 
                    top=0.90, 
                    wspace=1.0, 
                    hspace=0.0)
print ('Start Drawing' )
for iGrp in range(0,NumGrps):
    Y = [ 10.0 * Data[i][ilDrawGrpIdxs[iGrp] + 1] for i in range(NumData) ]

    axis[iGrp].plot(X, Y, color='blue',linewidth=0.5, label=slGroupsRaw[ilDrawGrpIdxs[iGrp]])
    axis[iGrp].set_title(slGroupsRaw[ilDrawGrpIdxs[iGrp]], fontsize=12, loc='center', y=0.60, bbox=dict(facecolor='white',alpha=0.65,pad=2,boxstyle='round,pad=0.1',linestyle='none'))
    axis[iGrp].set_xlim(xmin=t_ini,xmax=t_fin)

    yarr = np.asarray(Y)
    av[iGrp]  = np.mean(yarr)
    std[iGrp] = np.std(yarr)
    max[iGrp] = np.max(yarr)
    min[iGrp] = np.min(yarr)
    imax[iGrp] = yarr.argmax()
    imin[iGrp] = yarr.argmin()
    tmax = X[imax[iGrp]]
    tmin = X[imin[iGrp]]
    print ('%s: %4.3f +/- %5.4f, min = %4.3f, max = %4.3f , tmin=%7.3f, tmax=%7.3f' % ( slGroupsRaw[ilDrawGrpIdxs[iGrp]], av[iGrp], std[iGrp], min[iGrp], max[iGrp], tmin, tmax ) )
    print ('%s: %4.3f +/- %5.4f, min = %4.3f, max = %4.3f , tmin=%7.3f, tmax=%7.3f' % ( slGroupsRaw[ilDrawGrpIdxs[iGrp]], av[iGrp], std[iGrp], min[iGrp], max[iGrp], tmin, tmax ) , file = logFile)

    axis[iGrp].axhline(y=av[iGrp], color='b', linestyle='-', linewidth=0.7, label='mean value')

    if slGroupsRaw[ilDrawGrpIdxs[iGrp]] == 'REF:N20_REF:H30':
        continue
logFile.close()
print ('===== Done with all Groups =====' )


for ax in axis.flat:
    ax.label_outer()


fig_subpl.savefig('${dynfig}',dpi=200)
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




