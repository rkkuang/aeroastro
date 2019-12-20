from urllib import request
from urllib.request import urlopen
import json
import os
import sys
import time
import re
from subprocess import call
abspath = os.path.abspath('./')
savedir = "./v5_4_45"
headers ={
"User-Agent":"Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/77.0.3865.90 Safari/537.36"
}
baseurl = "http://www.sdss3.org/svn//repo/idlspec2d/tags/v5_4_45"

# course_url = "https://indico.skatelescope.org/event/551/overview"
# lecbase_url = "https://indico.skatelescope.org"
def get_url(course_url, headers):
    req = request.Request(url=course_url, headers = headers)#, method='POST'
    response = request.urlopen(req)
    content = response.read().decode('utf-8')
    # jcontent = json.loads(content)
    jcontent = content
    return jcontent

def dl(baseurl, savedir):
    jcontent = get_url(baseurl, headers)
    # pdflinks = re.findall("/event.*?.pdf", jcontent)
    dirlinks = re.findall("name=\"(.*?)\"", jcontent)
    dirlinks_full = re.findall("name=\"(.*)\"", jcontent)

    # print(dirlinks)
    # print(dirlinks_full)


    for l,lfull in zip(dirlinks,dirlinks_full):
        savedirtemp = savedir+"/"+l
        dirlink = baseurl+"/"+l
        print(dirlink)
        # dircontent = get_url(dirlink, headers)
        # if not "name=\"" in dircontent:
        if not "/" in lfull:
            if os.path.exists(savedirtemp) is False:
                cmd="wget {} -O {}".format(dirlink, savedirtemp)
                print(cmd)
                call(cmd.split())
        else:
            if os.path.exists(savedirtemp) is False:
                os.mkdir(savedirtemp)
            print("*"*25,dirlink, savedirtemp,"*"*25)
            dl(dirlink, savedirtemp)
            # filename = l.split("/")[-1]
            # filename = "{}/{}".format(savedir, filename)
            # if os.path.exists(savedir) is False:
            #     os.mkdir(savedir)
            # url = lecbase_url+l
            # # cmd="axel {} -o {}".format(url, filename)
            # cmd="wget {} -O {}".format(url, filename)
            # call(cmd.split())
dl(baseurl, savedir)