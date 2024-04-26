import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib.collections import LineCollection
from scipy.spatial.transform import Rotation as R
from scipy.optimize import minimize

# Импортируем класс слайдера
from matplotlib.widgets import Slider

def pyramid(a, b, h):
    '''А тут мы нарисуем пирамиду'''
    return np.array([[-a/2, -b/2, 0], [a/2, -b/2, 0], [a/2, b/2, 0],  [-a/2, b/2, 0], [0, 0, h]])

def doXZ(v):
    return v[:,[0,2]]

def rotate_pyramid(p, alpha, beta, gamma):
    r = R.from_euler('zyx', [alpha, beta, gamma], degrees=True)
    return r.apply(p)

def find_base_angle(p):
    return np.rad2deg(np.arccos(np.dot((p[1] - p[0])/np.linalg.norm(p[1] - p[0]), (p[3] - p[0])/np.linalg.norm(p[3] - p[0]))))
    
def find_2d_angle_from_rotation_minus_pi(params, visang):
    a, b, g = params
    pg = pyramid(2, 2, 2)
    pg = rotate_pyramid(pg, a, b, g)
    p = doXZ(pg)
    return np.abs(find_base_angle(p) - visang)

def find_orientation(projection):
    angl = find_base_angle(projection)
    res = minimize(find_2d_angle_from_rotation_minus_pi, [45, 1, 45], tol=1e-8, args=angl, method='Nelder-Mead', bounds=[(-90,90), (-1e-5, 1e-5), (-90, 90)])
    print(res.x)
    print(find_2d_angle_from_rotation_minus_pi(res.x, 0))
    
    

def updateGraph():
    '''!!! Функция для обновления графика'''
    global slider_a
    global slider_b
    global slider_h
    global slider_alpha
    global slider_beta
    global slider_gamma
    global graph_axes
    # Используем атрибут val, чтобы получить значение слайдеров
    a = slider_a.val
    b = slider_b.val
    h = slider_h.val
    # ВВращаем пирамиду
    v = rotate_pyramid(pyramid(a, b, h), slider_alpha.val, slider_beta.val, slider_gamma.val)
    # чистим оси
    graph_axes.clear()
    # рисуем пирамиду
    graph_axes.scatter3D(v[:, 0], v[:, 1], v[:, 2])
    verts = [ [v[0],v[1],v[4]], [v[0],v[3],v[4]], [v[2],v[1],v[4]], [v[2],v[3],v[4]], [v[0],v[1],v[2],v[3]]]
    graph_axes.add_collection3d(Poly3DCollection(verts, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.25))
    graph_axes.set_xticks([-6,-4,-2,0,2,4,6])
    graph_axes.set_yticks([-6,-4,-2,0,2,4,6])
    graph_axes.set_zticks([-6,-4,-2,0,2,4,6])
    graph_axes.axes.set_xlim3d(left=-6, right=6)
    graph_axes.axes.set_ylim3d(bottom=-6, top=6)
    graph_axes.axes.set_zlim3d(bottom=-6, top=6)
    
    # считаем проекцию на плоскость XZ
    p = doXZ(v)
    edges = [ [p[0], p[1]], [p[1], p[2]], [p[2], p[3]], [p[3], p[0]], [p[0], p[4]], [p[1], p[4]], [p[2], p[4]], [p[3], p[4]] ]
    lines = LineCollection(edges)
    proj_axes.clear()
    proj_axes.scatter(p[:,0], p[:,1])
    proj_axes.add_collection(lines)
    proj_axes.set_xlim([-6,6])
    proj_axes.set_ylim([-6,6])
    
    real_base_angl = find_base_angle(v)
    visible_base_angl = find_base_angle(p)
    print("real base angl = {}, visible base angl = {}".format(real_base_angl, visible_base_angl))
    find_orientation(p)

    plt.draw()


def onChangeValue(value: np.float64):
    '''!!! Обработчик события изменения значений слайдеров'''
    updateGraph()


if __name__ == '__main__':
    # Начальные параметры графиков
    current_sigma = 0.2
    current_mu = 0.0

    # Создадим окно с графиком
    fig = plt.figure()
    graph_axes = fig.add_subplot(221,projection='3d')
    proj_axes = fig.add_subplot(222, adjustable='box', aspect=1)
    #graph_axes.grid()

    # Выделим область, которую будет занимать график
    fig.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.4)

    axes_slider_a = plt.axes([0.05, 0.35, 0.75, 0.04])
    slider_a = Slider(axes_slider_a,
                          label='a',
                          valmin=0.1,
                          valmax=10.0,
                          valinit=2,
                          valfmt='%1.2f')

    axes_slider_b = plt.axes([0.05, 0.30, 0.75, 0.04])
    slider_b = Slider(axes_slider_b,
                       label='b',
                       valmin=0.1,
                       valmax=10.0,
                       valinit=2,
                       valfmt='%1.2f')

    axes_slider_h = plt.axes([0.05, 0.25, 0.75, 0.04])
    slider_h = Slider(axes_slider_h,
                       label='h',
                       valmin=0.1,
                       valmax=10.0,
                       valinit=2,
                       valfmt='%1.2f')

    axes_slider_alpha = plt.axes([0.05, 0.20, 0.75, 0.04])
    slider_alpha = Slider(axes_slider_alpha,
                          label='alpha',
                          valmin=-90,
                          valmax=90,
                          valinit=0,
                          valfmt='%1.2f')

    axes_slider_beta = plt.axes([0.05, 0.15, 0.75, 0.04])
    slider_beta = Slider(axes_slider_beta,
                       label='beta',
                       valmin=-90,
                       valmax=90,
                       valinit=0,
                       valfmt='%1.2f')

    axes_slider_gamma = plt.axes([0.05, 0.10, 0.75, 0.04])
    slider_gamma = Slider(axes_slider_gamma,
                       label='gamma',
                       valmin=-90,
                       valmax=90,
                       valinit=0,
                       valfmt='%1.2f')



    # !!! Подпишемся на события при изменении значения слайдеров.
    slider_a.on_changed(onChangeValue)
    slider_b.on_changed(onChangeValue)
    slider_h.on_changed(onChangeValue)

    slider_alpha.on_changed(onChangeValue)
    slider_beta.on_changed(onChangeValue)
    slider_gamma.on_changed(onChangeValue)


    updateGraph()
    plt.show()

