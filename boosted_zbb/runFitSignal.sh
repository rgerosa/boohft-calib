cd datacards  
text2workspace.py datacard_cat0.txt -m 90 -o datacard_cat0.root
text2workspace.py datacard_cat1.txt -m 90 -o datacard_cat1.root
text2workspace.py datacard_cat2.txt -m 90 -o datacard_cat2.root
text2workspace.py datacard_cat3.txt -m 90 -o datacard_cat3.root

text2workspace.py datacard_cat0_pass.txt -m 90 -o datacard_cat0_pass.root
text2workspace.py datacard_cat1_pass.txt -m 90 -o datacard_cat1_pass.root
text2workspace.py datacard_cat2_pass.txt -m 90 -o datacard_cat2_pass.root
text2workspace.py datacard_cat3_pass.txt -m 90 -o datacard_cat3_pass.root

text2workspace.py datacard_cat0_fail.txt -m 90 -o datacard_cat0_fail.root
text2workspace.py datacard_cat1_fail.txt -m 90 -o datacard_cat1_fail.root
text2workspace.py datacard_cat2_fail.txt -m 90 -o datacard_cat2_fail.root
text2workspace.py datacard_cat3_fail.txt -m 90 -o datacard_cat3_fail.root

cd -
combine -M FitDiagnostics -m 90 datacards/datacard_cat3_fail.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat3_fail
combine -M FitDiagnostics -m 90 datacards/datacard_cat3_pass.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat3_pass
combine -M FitDiagnostics -m 90 datacards/datacard_cat3.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat3

combine -M FitDiagnostics -m 90 datacards/datacard_cat2_fail.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat2_fail
combine -M FitDiagnostics -m 90 datacards/datacard_cat2_pass.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat2_pass
combine -M FitDiagnostics -m 90 datacards/datacard_cat2.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat2

combine -M FitDiagnostics -m 90 datacards/datacard_cat1_fail.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat1_fail
combine -M FitDiagnostics -m 90 datacards/datacard_cat1_pass.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat1_pass
combine -M FitDiagnostics -m 90 datacards/datacard_cat1.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat1

combine -M FitDiagnostics -m 90 datacards/datacard_cat0_fail.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat0_fail
combine -M FitDiagnostics -m 90 datacards/datacard_cat0_pass.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat0_pass
combine -M FitDiagnostics -m 90 datacards/datacard_cat0.root --setParameters r=1 --setParameterRanges r=[-10,10] --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveShapes --saveWithUncertainties -n _cat0

mkdir fit_results
mv higgsCombine*root fit_results
mv fitDiagnostics* fit_results
