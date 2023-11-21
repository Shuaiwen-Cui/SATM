# usage: in git bash, input sh update.sh - windows; or in terminal input ./update.sh and click approve git for windows pop up to run it
# usage: in terminal, input .\update.sh - mac / linux ? not sure

echo '--------upload files start--------'   
# enter the target folder
# cd ./

# git init
git add .
git status
git commit -m 'auto update SATM'
echo '--------commit successfully--------'

# git push -f https://github.com/Shuaiwen-Cui/SATM.git main
git push -u https://github.com/Shuaiwen-Cui/SATM.git main
# git remote add origin https://github.com/Shuaiwen-Cui/SATM.git
# git push -u origin main
echo '--------push to GitHub successfully--------'
