# -*- coding: utf-8 -*-
from __future__ import unicode_literals

#from django.shortcuts import render_to_response
from django.shortcuts import render
# Create your views here.
from django.template import RequestContext
from django.http import HttpResponse
from models import InputForm
from compute import SPDC ##compute
import os

def index(request):
    os.chdir(os.path.dirname(__file__))
    result = None
    if request.method == 'POST':   #POST
        form = InputForm(request.POST) #POST
        if form.is_valid():
            form2 = form.save(commit=False)
            result = SPDC(form2.wavep, form2.angle, form2.L, form2.distz, form2.wp)#(wavep, angle, L, distz, wp)
            result = result.replace('static/', '')
    else:
        form = InputForm()
   
   ## context ={'form':form,
   ##           'result': result,
   ##           }

    return render(request,'vibb.html',  #return render_to_response('vibb.html',
            {'form': form,
             'result': result,
             })      ## }, context_instance=RequestContext(request))
   ## return render(request, 'vibb.html',context)
