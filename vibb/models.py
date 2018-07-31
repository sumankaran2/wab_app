# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models

# Create your models here.
from django.forms import ModelForm
from math import pi

class Input(models.Model):
    wavep = models.FloatField(
        verbose_name=' Wavep (nm)', default=405.0)
    angle = models.FloatField(
        verbose_name='Thetap (deg.)', default=28.648)
    L = models.FloatField(
        verbose_name='L (mm)', default=2)
    distz = models.FloatField(
        verbose_name=' distz (cm)', default=100)
    wp = models.FloatField(
        verbose_name='Beam waist(um)', default=388)

class InputForm(ModelForm):
    class Meta:
        model = Input
        fields = "__all__"
