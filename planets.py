#!/usr/bin/env python3

# This file is part of Planets.
#
# Planets is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Planets is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Planets.  If not, see <http://www.gnu.org/licenses/>.

import pygame
import sys
from pygame.locals import *
import numpy as np
from random import random
from gravity import gravity

class Planets:
    # Length in AU; Time in days; Mass in Earth mass
    solarmass = 332946.04877
    earthmass = 5.97219e24 # kg
    km = 1. / 149597870.700
    s = 1. / 86400.

    def __init__(self):
        self.x = []
        self.v = []
        self.m = []
        self.r = []
        self.color = []
        self.image = []
        self.name = []
        self.rect = []
        self.circle = np.linspace(0., 2.*np.pi, 64)

    def add_planet(self, x=[0.,0.], v=None, m=1., r=4.e-5, e=0., color=None, image=None, name='', isun=0):
        '''
        x [AU]      : Position
        v [AU/day]  : Velocity
        m [M_Earth] : Planet mass
        r [AU]      : Planet radius
        e [1]       : Eccentricity of orbit
        color       : Draw color
        image       : Surface instance
        name        : Some name
        '''
        if len(self.x) == 0:
            self.x = np.asfortranarray(np.array(x).reshape((1,2)))
        else:
            self.x = np.asfortranarray(np.append(self.x, np.array(x).reshape((1,2)), axis=0))
        if v is None:
            v = self.Keplerian_orbit(e, m, isun=isun)
        if len(self.v) == 0:
            self.v = np.asfortranarray(np.array(v).reshape((1,2)))
        else:
            self.v = np.asfortranarray(np.append(self.v, np.array(v).reshape((1,2)), axis=0))
        self.m = np.append(self.m, m)
        self.r = np.append(self.r, r)
        self.image.append(image)
        if color is None:
            color = (255,255,255)
            if image is not None:
                color = image.get_at((int(image.get_width()/2), int(image.get_height()/2)))
        self.color.append(color)
        self.name.append(name)
        self.rect.append(None)

    def Keplerian_orbit(self, e, m, isun=0):
        aphelion = np.sum((self.x[-1,:] - self.x[isun,:])**2)**0.5
        a = aphelion / (1. + e)
        v = ((gravity.g * (self.m[isun] + m)) / a * (1. - e) / (1. + e))**0.5
        x, y = tuple(self.x[-1,:] - self.x[isun,:])
        alpha = np.arctan2(-y, x)
        return [v*np.sin(alpha), v*np.cos(alpha)]

    def step(self, dt, ndt):
        dt = gravity.integrate(ndt, dt, self.x, self.v, self.m, self.r, len(self.m))
        return dt

    def draw(self, screen, monitor_center, center, scale, show_orbits, tooltip):
        a, b, e, alpha, perihelion, aphelion, T, isun = gravity.orbit(self.x, self.v, self.m, len(self.m))
        isun -= 1

        for i in range(len(self.m)):
            if self.m[i] > 0:
                if 0 <= e[i] < 1.-1.e-6 and isun[i] >= 0 and (show_orbits or tooltip == i):
                    beta = alpha[i] + np.arctan2(*tuple(self.x[i,::-1] - self.x[isun[i],::-1]))
                    x, y = a[i] * np.cos(self.circle) + a[i] * e[i], \
                           b[i] * np.sin(self.circle)
                    x, y = x * np.cos(beta) - y * np.sin(beta), \
                           x * np.sin(beta) + y * np.cos(beta)
                    x, y = (self.x[isun[i],0] + x - center[0]) * scale + monitor_center[0], \
                           (self.x[isun[i],1] + y - center[1]) * scale + monitor_center[1]
                    pygame.draw.polygon(screen, self.color[i], list(zip(x,y)), 1)

        for i in range(len(self.m)):
            if self.m[i] > 0:
                x, y = tuple(np.int32((self.x[i,:]-center)*scale + monitor_center))
                rscaled = self.r[i] * scale * 1.e3
                rmin = 4
                if self.image[i] is None:
                    self.rect[i] = pygame.draw.circle(screen, self.color[i], (x, y), max(rmin, int(rscaled)))
                else:
                    dx, dy = self.image[i].get_width(), self.image[i].get_height()
                    r = 0.5 * max(dx, dy)
                    dx, dy = max(rmin*2, int(dx * rscaled/r)), max(rmin*2, int(dy * rscaled/r))
                    x, y = x - 0.5 * dx, y - 0.5 * dy
                    self.rect[i] = screen.blit(pygame.transform.scale(self.image[i], (dx,dy)), (x, y))

        if tooltip is not None:
            font = pygame.font.SysFont('Arial', 16, bold=False)
            message = [self.name[tooltip]]
            if perihelion[tooltip] > 0 or aphelion[tooltip] > 0:
                message.append('%.2f - %.2f AU'%(perihelion[tooltip], aphelion[tooltip]))
            if T[tooltip] > 0:
                message.append('T = %.2f days'%T[tooltip])
            for n,m in enumerate(message):
                text = font.render(m, True, self.color[tooltip])
                offset = np.array([5, -text.get_rect().height*(len(message)-n) - 5])
                screen.blit(text, tuple(np.int32((self.x[tooltip,:]-center)*scale + monitor_center + offset)))

class Universe:
    def __init__(self):
        pygame.init()
        monitor_height = pygame.display.Info().current_h
        monitor_width = pygame.display.Info().current_w
        self.screen = pygame.display.set_mode((monitor_width, monitor_height), pygame.HWSURFACE | pygame.DOUBLEBUF)
        self.screen.fill((0, 0, 0))
        self.surf = pygame.Surface((monitor_width, monitor_height))
        self.surf.fill((0,0,0))
        self.monitor_center = 0.5 * np.array([monitor_width, monitor_height])
        self.center = np.zeros(2)
        self.font = pygame.font.SysFont('Arial', 16, bold=False)

        pygame.display.set_caption('Planets')
        pygame.mouse.set_visible(True)
        pygame.key.set_repeat(200, 30)
        pygame.display.flip()
        self.clock = pygame.time.Clock()

        self.dt = None
        self.planets = Planets()
        self.solar_system()
        #self.neutron_star()

        self.zoom = min(100., 0.5*monitor_height / max(1.e-33, np.max(np.sum(self.planets.x**2, axis=1)**0.5)))
        self.zoom0 = self.zoom
        self.timestep = 1.
        self.timestep0 = self.timestep
        self.show_orbits = True

    def solar_system(self):
        self.planets.add_planet(x=[0.,0.], v=[0.,0.], m=Planets.solarmass, r=0.5*4.e4*Planets.km, color=(249,246,91), name='Sun')

        imgsurf = pygame.image.load('solar_system.png')

        planets = []
        # Semi-major axis [AU], eccentricity, mass [Earth masses], diameter [km]
        planets.append(( 0.387099, 0.205630,   0.055256,   4879., ( 13, 12,  0,  0), 'Mercury'))
        planets.append(( 0.723332, 0.006772,   0.815030,  12104., ( 31, 31,  0, 14), 'Venus'))
        planets.append(( 1.000002, 0.016709,   1.000000,  12756., ( 31, 31,  0, 47), 'Earth'))
        planets.append(( 1.523679, 0.093390,   0.107450,   6792., ( 16, 16,  0, 80), 'Mars'))
        planets.append(( 5.202600, 0.048500, 317.877470, 142984., (324,302,  0,105), 'Jupiter'))
        planets.append(( 9.554909, 0.055550,  95.162370, 120536., (415,454,344,  0), 'Saturn'))
        planets.append((19.218400, 0.046380,  14.534650,  51118., (117,118,767,  0), 'Uranus'))
        planets.append((30.110387, 0.009460,  17.145970,  49528., (113,113,771,124), 'Neptune'))

        for a,e,m,r,k,name in planets:
            x = a * (1. + e)
            alpha = 2.*np.pi * random()
            c = None
            img = None
            if len(k) == 3:
                c = k
            else:
                img = pygame.Surface(tuple(k[:2])).convert_alpha()
                img.fill((255,255,255,0))
                img.blit(imgsurf, tuple(-np.array(k[2:])))
            self.planets.add_planet(x=[self.planets.x[0,0] + x*np.cos(alpha), \
                                       self.planets.x[0,1] + x*np.sin(alpha)], m=m, e=e, r=0.5*r*Planets.km, color=c, image=img, name=name)

    def neutron_star(self):
        self.planets.add_planet(x=[-50.,2.], v=[300.*Planets.km/Planets.s, 0.], m=1.4*Planets.solarmass, r=12.*Planets.km, name='Neutron Star')

    def step(self):
        if self.dt is not None:
            ndt = int(round(self.timestep / self.dt))
        else:
            self.dt = 1.e-10
            ndt = 1
        time = self.dt * ndt
        self.dt = self.planets.step(self.dt, ndt)
        return time

    def draw(self, t, tooltip):
        self.planets.draw(self.screen, self.monitor_center, self.center, self.zoom, self.show_orbits, tooltip)
        x0, y0 = 20, 20
        text = self.font.render('%i days' % t, True, (255,255,255))
        self.screen.blit(text, (x0, y0))
        y0 += 30
        pygame.draw.line(self.screen, (255,255,255), (x0, y0), (x0+self.zoom, y0))

    def run(self):
        t = 0.
        while True:
            self.clock.tick(20)

            for event in pygame.event.get():
                if event.type == KEYDOWN:
                    if event.key == K_ESCAPE:
                        sys.exit(0)
                    elif event.key == K_PLUS:
                        self.timestep *= 2.
                    elif event.key == K_MINUS:
                        self.timestep *= 0.5
                    elif event.key == K_0:
                        self.timestep *= -1.
                    elif event.key == K_PAGEUP:
                        self.zoom *= 2.
                    elif event.key == K_PAGEDOWN:
                        self.zoom *= 0.5
                    elif event.key == K_HOME:
                        self.zoom = self.zoom0
                        self.center = np.zeros(2)
                        self.timestep = self.timestep0
                    elif event.key == K_LEFT:
                        self.center[0] = self.center[0] - 20 / self.zoom
                    elif event.key == K_RIGHT:
                        self.center[0] = self.center[0] + 20 / self.zoom
                    elif event.key == K_UP:
                        self.center[1] = self.center[1] - 20 / self.zoom
                    elif event.key == K_DOWN:
                        self.center[1] = self.center[1] + 20 / self.zoom
                    elif event.key == K_o:
                        self.show_orbits = not self.show_orbits

            self.screen.blit(self.surf, (0,0))

            tooltip = None
            for n,rect in enumerate(self.planets.rect):
                if rect is not None:
                    if rect.collidepoint(pygame.mouse.get_pos()):
                        tooltip = n

            if self.timestep > 0:
                t += self.step()
            self.draw(t, tooltip)
            pygame.display.flip()

u = Universe()
u.run()
