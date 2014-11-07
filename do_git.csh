set commitname = $argv[1]
git init
git config --global user.name "jiayixie"
git config --global user.email jiayi.seis@gmail.com
git add -A
git commit -m "$commitname"
git remote add origin git@github.com:jiayixie/inversion_ElasticTensor.git
git push -u origin master
