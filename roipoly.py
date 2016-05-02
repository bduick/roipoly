# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 21:15:44 2016

@author: duick
"""

# todo
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

# todo
# features to add:
# 1) in move mode if click segment, move the whole line
# 2) in add_verts mode, if click segment, add point on the segment
class roipoly:

    class mode:
        add_verts, move, finalize = range(3)

    def __init__(self, **kwargs):
        self.fig = plt.gcf()
        self.ax = plt.gca()
        # turn off auto-scale because adding plot on image changes image scale
        # even if plot is within the existing axes limits
        self.ax.set_autoscale_on(False)

        axes_images = self.ax.get_images()
        if len(axes_images) != 1:
            print("error")  # TODO error message
        self.image_h, self.image_w = axes_images[0].get_size()

        self.xs = [self.image_w/2]
        self.ys = [self.image_h/2]
        self.line, = self.ax.plot(self.xs, self.ys, '-ok', animated=True)
        self.line.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'),
                                    path_effects.Normal()])

        self.max_pixels_from_vertex = 5  # maximum pixel distance from vertex to consider as vertex selection
        self.x0t, self.y0t = None, None
        self.selected_vertex_idx = None
        self.mask = None
        self.current_mode = self.mode.add_verts

        canvas = self.fig.canvas
        # TODO read about draw_event and animated=True
        self.event_cids = []
        self.event_cids.append(canvas.mpl_connect('draw_event', self.draw_callback))
        self.event_cids.append(canvas.mpl_connect('button_press_event', self.button_press_callback))
        self.event_cids.append(canvas.mpl_connect('button_release_event', self.button_release_callback))
        self.event_cids.append(canvas.mpl_connect('motion_notify_event', self.motion_notify_callback))

        if 'block' in kwargs.keys() and kwargs['block'] == False:
            plt.show(block=False)
        else:
            plt.show()

    def draw_callback(self, event):
        # TODO figure out what these ax.bbox commands are
        self.background = self.fig.canvas.copy_from_bbox(self.ax.bbox)
        if self.line is not None:
            self.ax.draw_artist(self.line)
        self.fig.canvas.blit(self.ax.bbox)

    def button_press_callback(self, event):
        """TBD"""

        # do nothing if mouse event was not in the image axes
        if event.inaxes != self.ax:
            return

        # switch modes if right-mouse click
        if event.button == 3:
            if self.current_mode == self.mode.add_verts and len(self.xs) > 1:
                self.current_mode = self.mode.move
                self.xs.pop()
                self.ys.pop()
                self.line.set_data(self.xs, self.ys)
                self.line.set_color('red')
            elif self.current_mode == self.mode.move:
                if len(self.xs) > 2:
                    self.current_mode = self.mode.finalize
                    self.xs.append(self.xs[0])
                    self.ys.append(self.ys[0])
                    self.line.set_data(self.xs, self.ys)
                    self.line.set_color('blue')
                else:
                    self.xs.append(event.xdata)
                    self.ys.append(event.ydata)
                    self.line.set_data(self.xs, self.ys)
                    self.line.set_color('black')
                    self.current_mode = self.mode.add_verts
            elif self.current_mode == self.mode.finalize:
                self.xs[-1] = event.xdata
                self.ys[-1] = event.ydata
                self.line.set_data(self.xs, self.ys)
                self.line.set_color('black')
                self.current_mode = self.mode.add_verts
            # TODO consider moving line.set_data down here; also mutex lock
            self.fig.canvas.draw()
            return

        # do nothing if not left-mouse click
        if event.button != 1:
            return

        # add vertex mode
        if self.current_mode == self.mode.add_verts:
            if len(self.xs) > 3 and self.xs[0] == self.xs[-1] and self.ys[0] == self.ys[-1]:
                self.current_mode = self.mode.finalize
                return
            self.xs.insert(-1, event.xdata)
            self.ys.insert(-1, event.ydata)
            # TODO add mutex lock
            self.line.set_data(self.xs, self.ys)
            self.fig.canvas.draw()

            # TODO document what this is for
            self.x0t, self.y0t = self.line.get_transform().transform((self.xs[0], self.ys[0]))

        # move poly mode
        if self.current_mode == self.mode.move:
            self.selected_vertex_idx = self.get_index_under_point(event)

        # calculate ROI mask and exit if double click TODO wording
        if self.current_mode == self.mode.finalize and event.dblclick and \
                self.line.get_path().contains_point((event.xdata, event.ydata)):
            self.calc_mask()
            self.exit()
            return

    def button_release_callback(self, event):
        self.selected_vertex_idx = None
        # TODO document what this is for
        self.x0t, self.y0t = self.line.get_transform().transform((self.xs[0], self.ys[0]))

    def motion_notify_callback(self, event):
        """on mouse movement"""
        if event.inaxes is None:
            return

        # TODO combine these two functions

        if self.current_mode == self.mode.add_verts:
            self.xs[-1], self.ys[-1] = event.xdata, event.ydata
            self.line.set_data(self.xs, self.ys)

            if event.button is None and len(self.xs) > 3:
                d = np.sqrt((self.x0t - event.x)**2 + (self.y0t - event.y)**2)
                if d < self.max_pixels_from_vertex:
                    self.line.set_color('b')
                    self.xs[-1], self.ys[-1] = self.xs[0], self.ys[0]
                else:
                    self.line.set_color('k')

            self.fig.canvas.restore_region(self.background)
            self.ax.draw_artist(self.line)
            self.fig.canvas.blit(self.ax.bbox)
            return

        if self.current_mode == self.mode.move and self.selected_vertex_idx is not None:
            idx = self.selected_vertex_idx
            self.xs[idx], self.ys[idx] = event.xdata, event.ydata
            self.line.set_data(self.xs, self.ys)
            self.fig.canvas.restore_region(self.background)
            self.ax.draw_artist(self.line)
            self.fig.canvas.blit(self.ax.bbox)
            return

    def get_index_under_point(self, event):
        """get the index of the vertex under point if within epsilon tolerance"""
        xy = np.asarray(list(zip(self.xs, self.ys)))
        xyt = self.line.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.sqrt((xt - event.x) ** 2 + (yt - event.y) ** 2)
        pt_idx = np.argmin(d)
        if d[pt_idx] >= self.max_pixels_from_vertex:
            pt_idx = None
        return pt_idx

    def calc_mask(self):
        image_coord_x = np.tile(np.arange(self.image_w), self.image_h)
        image_coord_y = np.repeat(np.arange(self.image_h), self.image_w)
        image_coords = np.vstack((image_coord_x,image_coord_y)).T
        # TODO clean up
        result = self.line.get_path().contains_points(image_coords)
        self.mask = result.reshape((self.image_h, self.image_w))

    def exit(self):
        plt.close(self.fig)


if __name__ == '__main__':
    # sample_img = np.tile(np.arange(0, 1, 0.1), [15, 1])
    sample_img = np.load('sample_img.npy')
    plt.imshow(sample_img)
    roi = roipoly()
    mask = roi.mask
    if mask is not None:
        mask3 = np.tile(np.expand_dims(mask,axis=2), [1, 1, 3])
        plt.imshow(mask3*sample_img)
        plt.show()