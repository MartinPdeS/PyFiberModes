from PyFiberModes import FiberFactory, HE11, HE12, LP01, LP02, Wavelength, Mode
from PyFiberModes.field import Field
from MPSPlots.render2D import SceneList

# if __name__ == '__main__':
#     factory = FiberFactory()

#     factory.addLayer(
#         name='core',
#         radius=4.5e-6,
#         index=1.4489,
#         geometry="StepIndex"
#     )

#     factory.addLayer(
#         name='clad',
#         radius=62.5e-6,
#         index=1.4440,
#         geometry="StepIndex"
#     )

#     # factory.addLayer(
#     #     name='air',
#     #     radius=1e-6,
#     #     index=1.0001
#     # )

#     fiber = factory[0]

#     # print(factory._fibers)
#     # dsa

#     figure = SceneList(unit_size=(4, 4))

#     for mode in [HE11, LP01]:
#         field = Field(
#             fiber=fiber,
#             mode=mode,
#             wavelength=1550e-9,
#             limit=100e-6
#         )

#         ex = field.Ex()

#         ax = figure.append_ax(title=f'{mode}')

#         ax.add_mesh(
#             scalar=ex
#         )

#     figure.show()
