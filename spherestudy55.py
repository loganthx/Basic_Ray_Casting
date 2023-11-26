import numpy as np, matplotlib.pyplot as plt
from tqdm import tqdm

# w, h = 1024, 768
w, h = 160, 120
canvas = np.zeros((h, w, 3))
# you paint the canvas like canvas[y][x] 


class Sphere:
	def __init__(self, radius, center, color, specular):
		self.radius = radius
		self.center = center
		self.color = color
		self.specular = specular


class PointLight:
	def __init__(self, origin, intensity):
		self.origin = origin 
		self.intensity = intensity
		self.type = 'point'

class AmbientLight:
	def __init__(self, intensity):
		self.intensity = intensity	
		self.type = 'ambient'

class DirectionalLight:
	def __init__(self, direction, intensity):
		self.direction = direction
		self.intensity = intensity		
		self.type = 'directional'


def compute_diffuse_light_intensity(P, N):
	i = 0.0
	for light in scene_lights:
		if light.type == 'ambient':
			i += light.intensity

		if light.type == 'point':
			L_point = light.origin - P
			n_dot_l = np.dot(N, L_point)
			if n_dot_l > 0:
				length_N = np.sqrt(N[0]**2 + N[1]**2 + N[2]**2)
				lenght_L = np.sqrt( L_point[0]**2 + L_point[1]**2 + L_point[2]**2)
				i += light.intensity * n_dot_l / (length_N * lenght_L)

		if light.type == 'directional':
			L_directional = light.direction
			n_dot_l = np.dot(N, L_directional)
			if n_dot_l > 0:
				length_N = np.sqrt(N[0]**2 + N[1]**2 + N[2]**2)
				lenght_L = np.sqrt(L_directional[0]**2 + L_directional[1]**2+ L_directional[2]**2)
				i += light.intensity * n_dot_l / (length_N * lenght_L)

	return i




def compute_specular_light_intensity(sphere, P, N):
	i = 0.0
	# Define Eye Vector BeforeHand:
	V = -D; length_V = np.sqrt(V[0]**2 + V[1]**2 + V[2]**2)

	for light in scene_lights:

		if light.type == 'point':
			L_point = light.origin - P
			R = 2 * N * np.dot(N, L_point) - L_point
			r_dot_v = np.dot(R, V)
			if r_dot_v > 0:
				length_R = np.sqrt(R[0]**2 + R[1]**2 + R[2]**2)
				i += light.intensity * np.power(r_dot_v / (length_R*length_V), sphere.specular)

		if light.type == 'directional':
			L_directional = light.direction
			R = 2 * N * np.dot(N, L_directional) - L_directional
			r_dot_v = np.dot(R, V)
			if r_dot_v > 0:
				length_R = np.sqrt(R[0]**2 + R[1]**2 + R[2]**2)
				i += light.intensity * np.power(r_dot_v/ (length_R*length_V), sphere.specular)

	return i




def find_t(O, D, sphere):
	r = sphere.radius; CO = O - sphere.center
	a = np.dot(D, D); b = 2*np.dot(CO, D); c = np.dot(CO, CO) - r*r
	discriminant = b*b - 4*a*c
	if discriminant >= 0:
		t1 = (-b + np.sqrt(discriminant))/(2*a)
		t2 = (-b - np.sqrt(discriminant))/(2*a)
	else:
		t1, t2 = False, False
	return t1, t2


def paint_canvas(P, sphere, light):
	x = P[0]*w/P[2] + w/2
	y = h/2 - P[1]*h/P[2]
	x, y = int(x), int(y)
	color = sphere.color 
	try: 
		canvas[y][x] = light*color
	except: pass

O = np.array([0, 0, 0])
sphere_one = Sphere(radius = 1, center = np.array([0, -1, 3]), color = np.array([0, 0, 1]), specular = 500)
sphere_two = Sphere(radius = 1, center = np.array([2, 0, 4]), color = np.array([1, 0, 1]), specular = 50)
sphere_three = Sphere(radius = 1, center = np.array([-2, 0, 4]), color = np.array([1, 0, 0]), specular = 100)
sphere_four = Sphere(radius = 5000, center = np.array([0, -5001, 0]), color = np.array([1, 1, 0]), specular = 1000)

point_light = PointLight(origin = np.array((0, 0, 0)), intensity = .5)
ambient_light = AmbientLight(intensity = .3)
directional_light = DirectionalLight(direction = np.array((1, 4, 4)), intensity = .2)

scene_lights = [point_light, ambient_light, directional_light]
scene_spheres = [sphere_four, sphere_one, sphere_two, sphere_three]


for sphere in tqdm(scene_spheres):
	for x in np.arange(-w/2, w/2, .5):
		for y in np.arange(-h/2, h/2, .5):
			D = np.array([x/w, y/h, 1])
			
			t1, t2 = find_t(O, D, sphere)
			if min(t1, t2):
				t = min(t1, t2)
				if t > 1:
					P = O + t*D
					N = P - sphere.center 
					diffuse = compute_diffuse_light_intensity(P, N)
					specular = compute_specular_light_intensity(sphere, P, N)
					light = np.clip( (diffuse + specular ),  0, 1 )
					paint_canvas(P, sphere, light)


plt.imshow(canvas)
plt.show()