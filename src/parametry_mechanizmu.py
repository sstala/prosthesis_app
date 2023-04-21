# %matplotlib qt

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plotly.graph_objects as go
import plotly
import plotly_express as px
import plotly.io as pio
import plotly.express as px
pio.renderers.default='notebook_connected'
# pio.renderers.default='browser'



def auxliary_edge_length(p1x,p2x,p1y,p2y):
    return ((p1x-p2x)**2+(p1y-p2y)**2)**0.5
    
def triangle_perimeter(len1,len2,len3):
    return (len1+len2+len3)/2

def triangle_area(perimeter,len1,len2,len3):
    return (perimeter*(perimeter-len1)*(perimeter-len2)*(perimeter-len3))**0.5

# def find_dimensions():
#     mx = 0
#
#     stp = 0.01
#     nx = np.arange(-10,10+stp,stp)
#     ny = np.arange(-8,8+stp,stp)
#     l1 = 35
#
#     my = np.zeros((len(nx), len(ny)))
#     for i in range(0,len(nx)):
#         for j in range(0,len(ny)):
#             my[i,j] = (-l1*(mx+nx[i]+ny[j]) + 2*mx*nx[i])/(-l1-2*ny[j])
#
# fig = go.Figure(go.Surface(
#     x = nx,
#     y = ny,
#     z = my
#     ))
# fig.update_layout(
#         scene = {
#             "xaxis": {"nticks": 20},
#             "zaxis": {"nticks": 4},
#             'camera_eye': {"x": 0, "y": -1, "z": 0.5},
#             "aspectratio": {"x": 1, "y": 1, "z": 0.2}
#         })
#
# fig.show()
# Obliczanie parametrów mechanizmu

Ax = 0
Ay = 0
AB = 40
Ex = 0
Ey = -5.5

Fx0 = 40 - 5.5
Fy0 = 0
Bx0 = 40
By0 = 0
Gx0 = 40
Gy0 = -5.5
Hx0 = 40 + 35 - 3.5
Hy0 = 3.1
Cx0 = 40 + 35
Cy0 = 0
Dx0 = 40 + 35 + 25
Dy0 = 0

BF = 5.5
EF = (34.5**2+5.5**2)**0.5
BC = 35


AE = auxliary_edge_length(Ax,Ex,Ay,Ey)


joint_trajectory = pd.DataFrame(columns=['Ax', 'Ay','Bx', 'By','Cx', 'Cy','Dx', 'Dy','Ex', 'Ey','Fx', 'Fy','Gx', 'Gy','Hx', 'Hy'])
# mechanism_position
for i in np.linspace(0,90,91):
    alfa = i
    Bx = Ax + AB*np.cos(alfa*np.pi/180)
    By = Ay - AB*np.sin(alfa*np.pi/180)
    
    BE = auxliary_edge_length(Ex,Bx,Ey,By)
    
    per0 = triangle_perimeter(AE,AB,BE)
    are0 = triangle_area(per0,AE,AB,BE)
    gam0 = 90 + alfa - np.arcsin(are0/(0.5*AB*BE)) * (180/np.pi)
    
    per1 = triangle_perimeter(BF,EF,BE)
    are1 = triangle_area(per1,BF,EF,BE)
    gam1 = np.arcsin(are1/(0.5*EF*BE)) * (180/np.pi)
    
    Fx = Ex + EF * np.sin((gam0-gam1)*np.pi/180)
    Fy = Ey + EF * np.cos((gam0-gam1)*np.pi/180)
    
    tet1 = np.arccos(np.round_(((Bx-Fx)/BF),4)) * 180/np.pi
    
    BG = 5.5
    
    Gx = Bx - BG*np.sin(alfa*np.pi/180)
    Gy = By - BG*np.cos(alfa*np.pi/180)
    
    Cx = Bx + BC * np.cos(tet1*np.pi/180)
    Cy = By - BC * np.cos((90-tet1)*np.pi/180)
    
    tet2 = np.arccos(np.round_(((Cx-Bx)/BC),4)) * 180/np.pi
    
    CG = auxliary_edge_length(Cx,Gx,Cy,Gy)
    per2 = triangle_perimeter(BG,BC,CG)
    are2 = triangle_area(per2,BG,BC,CG)
    gam21 = np.arcsin(are2/(0.5*BC*CG)) * 180/np.pi
    gam22 = np.arcsin(are2/(0.5*BG*BC)) * 180/np.pi
    gam2 = 180 - gam21 - gam22
    
    GH = ((5.5 + 3.5)**2 + (35-3.1)**2)**0.5
    CH = (3.1**2 + 3.5**2)**0.5
    
    per3 = triangle_perimeter(GH,CH,CG)
    are3 = triangle_area(per3,GH,CH,CG)
    gam3 = np.arcsin(are3/(0.5*GH*CG)) * 180/np.pi
    
    tet3 = gam2 - gam3 + alfa
    
    Hx = Gx + GH * np.sin(tet3*np.pi/180)
    Hy = Gy + GH * np.cos(tet3*np.pi/180)
    
    BH = auxliary_edge_length(Bx, Hx, By, Hy)
    
    per4 = triangle_perimeter(BH, CH, BC)
    are4 = triangle_area(per4, BH, CH, BC)
    gam4 = np.arccos((CH**2 + BC**2 - BH**2) / (2 * CH * BC)) * 180/np.pi
    
    CD = 25
    DH = ((25+3.1)**2 + 3.5**2)**0.5
    
    per5 = triangle_perimeter(CH, CD, DH)
    are5 = triangle_area(per5, CH, CD, DH)
    gam5 = np.arccos((CH**2 + CD**2 - DH**2) / (2 * CH * CD)) * 180/np.pi
    
    tet4 = gam5 + gam4 - (90 - tet1)
    
    Dx = Cx + CD * np.sin(tet4 * np.pi/180)
    Dy = Cy + CD * np.cos(tet4 * np.pi/180)

    new_row = {'Ax': Ax, 'Ay': Ay,'Bx': Bx, 'By': By,'Cx': Cx, 'Cy': Cy,'Dx': Dx, 'Dy': Dy,'Ex': Ex, 'Ey': Ey,'Fx': Fx, 'Fy': Fy,'Gx': Gx, 'Gy': Gy,'Hx': Hx, 'Hy': Hy}
    joint_trajectory = joint_trajectory.append(new_row, ignore_index=True)
    
    


plt.figure(figsize = (12,6))

plt.plot(joint_trajectory['Bx'],joint_trajectory['By']) # trajektoria przegubów
plt.plot(joint_trajectory['Cx'],joint_trajectory['Cy']) # trajektoria przegubów
plt.plot(joint_trajectory['Dx'],joint_trajectory['Dy']) # trajektoria przegubów

plt.plot(joint_trajectory[['Ax', 'Bx', 'Cx', 'Dx']].iloc[0,:],joint_trajectory[['Ay', 'By', 'Cy', 'Dy']].iloc[0,:])
plt.plot(joint_trajectory[['Ax', 'Bx', 'Cx', 'Dx']].iloc[45,:],joint_trajectory[['Ay', 'By', 'Cy', 'Dy']].iloc[45,:])
plt.plot(joint_trajectory[['Ax', 'Bx', 'Cx', 'Dx']].iloc[90,:],joint_trajectory[['Ay', 'By', 'Cy', 'Dy']].iloc[90,:])
plt.xlim([-60,120])
plt.ylim([-80, 10])
plt.grid()
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
plt.show()

# aaa = joint_trajectory[['Ax', 'Bx', 'Cx', 'Dx']].iloc[1,:]


# import dash
# import dash_core_components as dcc
# import dash_html_components as html
# from dash.dependencies import Input, Output
#
# app = dash.Dash()
#
# app.layout = html.Div([
#     html.Div([
#         html.H2("Parametry mechanizmu")
#     ], className="banner"),
#
#     html.Div([
#         # html.Div([
#         #     html.H3("Wykres"),
#         #     dcc.Graph(
#         #         id = 'id01',
#         #         figure = {
#         #             'data': [
#         #                 {'x': joint_trajectory[['Ax', 'Bx', 'Cx', 'Dx']].iloc[0,:], 'y': joint_trajectory[['Ay', 'By', 'Cy', 'Dy']].iloc[0,:], 'type': 'line', 'line': {'color': '#E03A53', 'width': 5}, 'marker': {'color': '#E0AD0D'}},
#         #                 {'x': joint_trajectory[['Ax', 'Bx', 'Cx', 'Dx']].iloc[45,:], 'y': joint_trajectory[['Ay', 'By', 'Cy', 'Dy']].iloc[45,:], 'type': 'line', 'line': {'color': '#E03A53'}, 'marker': {'color': '#E0AD0D'}},
#         #                 {'x': joint_trajectory[['Ax', 'Bx', 'Cx', 'Dx']].iloc[90,:], 'y': joint_trajectory[['Ay', 'By', 'Cy', 'Dy']].iloc[90,:], 'type': 'line', 'line': {'color': '#E03A53'}, 'marker': {'color': '#E0AD0D'}},
#
#         #             ],
#         #             'layout': {
#         #                 'grid': {
#         #                     'color': '#ffffff'
#         #                 },
#         #                 'griddash': 'dot',
#         #                 'gridcolor': '#ffffff' ,
#         #                 'gridwidth': 3,
#         #                 'font': {
#         #                     'color': '#ffffff'
#         #                 },
#         #                 'plot_bgcolor': '#1e2130',
#         #                 'paper_bgcolor': '#1e2130',
#         #                 'margin': {
#         #                     'l': 30, 'r': 30, 't': 30, 'b': 30
#         #                 },
#         #                 'showlegend': False,
#         #                 'width': 540, 'height': 300,
#         #                 'xaxis': {
#         #                     'range': list([-20,100]),
#         #                     'gridcolor': '#70735d' ,
#         #                     'linecolor': '#70735d' ,
#         #                     'linewidth': 1,
#         #                     'mirror': True,
#         #                     'gridwidth': 0.05,
#         #                     'tickmode': 'linear', 'tick0': -60, 'dtick': 20,
#         #                     'scaleanchor': "y",
#         #                     'scaleratio': 1,
#
#         #                 },
#         #                 'yaxis': {
#         #                     'range': list([-80,20]),
#         #                     'scaleanchor': "x",
#         #                     'scaleratio': 1,
#         #                     'gridcolor': '#70735d' ,
#         #                     'linecolor': '#70735d' ,
#         #                     'linewidth': 1,
#         #                     'mirror': True,
#         #                     'gridwidth': 0.05,
#         #                     'tickmode': 'linear', 'tick0': 20, 'dtick': 20,
#         #                     'zeroline': False,
#
#         #                 },
#
#         #             }
#         #         }
#         #     )
#         # ], className="six columns"),
#         # html.Div([
#         #     html.H3("Inne")
#         # ], className="six columns")
#
#         html.Div(
#             id = "tabs",
#             className = "tabs",
#             children = [
#                 dcc.Tabs(
#                     id = "app-tabs",
#                     value = "tab_default",
#                     className = "custom-tabs",
#                     children = [
#                         dcc.Tab(
#                             id = "tab0",
#                             label = "INFO",
#                             value = "tab0",
#                             className = "custom-tab",
#                             selected_className = "custom-tab--selected",
#                         ),
#                         dcc.Tab(
#                             id = "tab1",
#                             label = "ZNAJDŹ WYMIARY",
#                             value = "tab1",
#                             className = "custom-tab",
#                             selected_className = "custom-tab--selected",
#                         ),
#                         dcc.Tab(
#                             id = "tab2",
#                             label = "POKAŻ NA WYKRESIE",
#                             value = "tab2",
#                             className = "custom-tab",
#                             selected_className = "custom-tab--selected",
#                         ),
#                         dcc.Tab(
#                             id = "tab3",
#                             label = "inne",
#                             value = "tab3",
#                             className = "custom-tab",
#                             selected_className = "custom-tab--selected",
#                         ),
#                     ]
#                 )
#
#             ]
#
#
#         ),
#
#         # html.Div(
#         #     id = 'container',
#         #     className = "big-container"
#         # )
#     ])
# ])
#
# # def show_info():
# #     return html.H3("Informacje")
#
# # def show_find_dimensions():
# #     # find_dimensions()
# #     return html.H3("znajdz wymiary")
#
# # def show_in_chart():
# #     return html.H3("poka sowe")
#
# # def show_other():
# #     return html.H3("inne")
#
#
#
# # @app.callback(
# #     Output("container", "children"),
# #     Input("app-tabs", "value")
# # )
#
# # def render_content(tab):
#
# #     if tab == "tab0":
# #         return show_info()
# #     elif tab == "tab1":
# #         return show_find_dimensions()
# #     elif tab == "tab2":
# #         return show_in_chart()
# #     elif tab == "tab3":
# #         return show_other()
#
#
# if __name__ == "__main__":
#     app.run_server(port=4050)