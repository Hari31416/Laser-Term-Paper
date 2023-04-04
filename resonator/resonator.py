import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.size"] = 16


class Mirror:
    """A class for a mirror object"""

    THETA1 = 145
    THETA2 = 215
    MIRROR_TICKS_COLOR = "k"

    def __init__(self, R, x, y, color=None, name="") -> None:
        """Instantiates the mirror object

        Parameters
        ----------
        R : float
            The radius of the mirror
        x : float
            The x-coordinate of the mirror
        y : float
            The y-coordinate of the mirror
        color : str, optional
            The color of the mirror, by default None
        name : str, optional
            The name of the mirror, by default ""
        """
        self.R = -R
        self.x = x
        self.y = y
        self.name = name
        self.theta1 = self.THETA1
        self.theta2 = self.THETA2
        self.color = color
        self.__setup()

    def __setup(self):
        """Sets up the mirror object"""
        mirror_up_point = np.abs(self.R) * np.cos(np.deg2rad(self.theta1))
        mirror_down_point = np.abs(self.R) * np.cos(np.deg2rad(self.theta2))
        mirror_length = mirror_up_point + mirror_down_point
        self.mirror_length = abs(mirror_length / 2)

    def get_x(self, y):
        """Returns the x-coordinate of the mirror at a given y-coordinate

        Parameters
        ----------
        y : float
            The y-coordinate

        Returns
        -------
        float
        """
        R = -self.R
        if self.x < 0:
            sintheta = (y - self.y) / R
            x = self.x + R * np.cos(np.arcsin(sintheta))
            return x - R
        else:
            sintheta = (y - self.y) / R
            x = self.x - R * np.cos(np.arcsin(sintheta))
            return x + R

    def plot(self, ax, hline=True):
        """Plots the mirror along with the mirror ticks

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes object
        hline : bool, optional
            Whether to plot the horizontal line passing through the mirror, by default True

        Returns
        -------
        None
        """
        theta1 = self.theta1
        theta2 = self.theta2
        thetas = np.linspace(theta1, theta2, 100)
        a, b = self.R, self.R

        if self.x < 0:
            x = self.x + a
            xs = a * np.cos(np.deg2rad(thetas)) + x
            ys = b * np.sin(np.deg2rad(thetas)) + self.y

        else:
            x = self.x - a
            xs = -a * np.cos(np.deg2rad(thetas)) + x
            ys = -b * np.sin(np.deg2rad(thetas)) + self.y

        ax.plot(xs, ys, color=self.color, linewidth=2)
        if hline:
            ax.axvline(x=self.x, color="k", linestyle="--")

        angles = np.arctan2(ys, xs)
        angles = np.abs(angles)
        if a < 0:
            mask = ys >= self.y

        else:
            mask = ys <= self.y
        angles[mask] = -angles[mask]

        for i in range(0, len(xs), 4):
            p1 = [xs[i], xs[i] + 0.1 * np.cos((angles[i]))]
            p2 = [ys[i], ys[i] + 0.1 * np.sin((angles[i]))]
            ax.plot(p1, p2, color=self.MIRROR_TICKS_COLOR, linewidth=1)


class Resonator:
    """A class to depict the resonator object"""

    LEFT_MIRROR_COLOR = "b"
    RIGHT_MIRROR_COLOR = "g"
    RAY_COLOR = "r"

    def __init__(self, R1, R2, L, name="") -> None:
        """Instantiates the resonator object

        Parameters
        ----------
        R1 : float
            The radius of the first circle
        R2 : float
            The radius of the second circle
        L : float
            The length of the resonator
        name : str, optional
            The name of the resonator, by default ""
        """
        self.R1 = R1
        self.R2 = R2
        self.L = L
        self.name = name
        self.g1 = 1 + self.L / self.R1
        self.g2 = 1 + self.L / self.R2
        self.fig = None
        self.ax = None
        self.y0 = None
        self.x0 = None
        self.__setup()

    def __setup(self):
        """Sets up the resonator object"""
        y_poistion = 0
        mirror_one_x = -self.L / 2
        mirror_two_x = self.L / 2
        self.left_mirror = Mirror(
            self.R1, mirror_one_x, y_poistion, color=self.LEFT_MIRROR_COLOR
        )
        self.right_mirror = Mirror(
            self.R2, mirror_two_x, y_poistion, color=self.RIGHT_MIRROR_COLOR
        )
        self.translation_rtm = np.array([[1, self.L], [0, 1]])
        self.left_reflection_rtm = np.array([[1, 0], [2 / self.R1, 1]])
        self.right_reflection_rtm = np.array([[1, 0], [2 / self.R2, 1]])

    @property
    def is_stable(self):
        """Checks if the resonator is stable"""
        if self.g1 * self.g2 > 1:
            return False
        if self.g1 * self.g2 < 0:
            return False
        return True

    def __metadata(self):
        """Returns the metadata of the resonator object"""
        is_stable = self.is_stable
        if is_stable:
            is_stable = "Stable"
        else:
            is_stable = "Unstable"
        return {
            "R1": round(self.R1, 2),
            "R2": round(self.R2, 2),
            "L": round(self.L, 2),
            "g1": round(self.g1, 2),
            "g2": round(self.g2, 2),
            "stability": is_stable,
        }

    def _metadata(self):
        """Returns the metadata of the resonator object as a string"""
        metadata = self.__metadata()
        metadata_text = ""
        i = 0
        for key, value in metadata.items():
            if i in [0, 1, 2]:
                metadata_text += f"{key} : {value:.2f} m\n"
            elif i in [3, 4]:
                metadata_text += f"{key} : {value:.2f}\n"
            else:
                metadata_text += f"{key} : {value}"
            i += 1
        return metadata_text

    def draw(
        self,
        zero_line=True,
        show=True,
        show_metadata=True,
    ):
        """Draws the resonator

        Parameters
        ----------
        zero_line : bool, optional
            Whether to show the zero line, by default True
        show : bool, optional
            Whether to show the plot, by default True
        show_metadata : bool, optional
            Whether to show the metadata, by default True

        Returns
        -------
        None
        """
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        self.left_mirror.plot(ax)
        self.right_mirror.plot(ax)

        offset_x_left = min(self.L / 2 * 0.8, 0.6)
        offset_x_right = min(self.L / 2 * 0.8, 0.6)
        if show_metadata:
            offset_x_right += 0.35
        R = max(np.abs(self.R1), np.abs(self.R2))
        x_ticks = np.round(
            np.linspace(-self.L / 2 - offset_x_left, self.L / 2 + offset_x_right, 11), 2
        )
        y_ticks = np.round(np.linspace(-R * 0.8, R * 0.8, 11), 2)

        ax.set_xlim(-self.L / 2 - offset_x_left, self.L / 2 + offset_x_right)
        ax.set_ylim(-R * 0.8, R * 0.8)
        ax.set_yticks(y_ticks)
        ax.set_xticks(x_ticks)
        ax.set_xlabel("x in m")
        ax.set_ylabel("y in m")

        if show_metadata:
            metadata = self._metadata()
            ax.text(
                0.72,
                0.75,
                metadata,
                transform=ax.transAxes,
                verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
            )

        if zero_line:
            ax.axhline(y=0, color="k", linestyle="--", alpha=0.5)

        ax.set_title(self.name)
        if show:
            plt.show()

        self.fig = fig
        self.ax = ax

    def get_left_x(self, y):
        """Returns the x position of the ray on the left mirror\\
        Uses the left mirror's get_x method"""
        return self.left_mirror.get_x(y)

    def get_right_x(self, y):
        """Returns the x position of the ray on the right mirror\\
        Uses the right mirror's get_x method"""
        return self.right_mirror.get_x(y)

    def __propogate_one(self, pos):
        """Propogates the ray through the resonator

        Parameters
        ----------
        pos : np.array
            The position vector of the ray

        Returns
        -------
        (yl, yr), pos5 : tuple
            The y position of the ray on the left and right mirrors and the final position vector
        """
        pos2 = self.translation_rtm @ pos
        yr = pos2[0]
        pos3 = self.right_reflection_rtm @ pos2
        assert np.isclose(pos3[0], yr), "Left y does not match"
        pos4 = self.translation_rtm @ pos3
        pos5 = self.left_reflection_rtm @ pos4
        yl = pos5[0]
        assert np.isclose(pos5[0], yl), "Right y does not match"
        return (yl, yr), pos5

    def __propogate(self, pos0, n):
        """Propogates the ray through the resonator

        Parameters
        ----------
        pos0 : np.array
            The position vector of the ray
        n : int, optional
            The number of times to propogate the ray

        Returns
        -------
        yls, yrs : tuple
            The y position of the ray on the left and right mirrors
        """
        yls = [pos0[0]]
        yrs = []
        for i in range(n):
            (yl, yr), pos0 = self.__propogate_one(pos0)
            yls.append(yl)
            yrs.append(yr)
        return yls, yrs

    def __preprocess(self, yls, yrs, n):
        """Preprocesses the y positions of the ray on the mirrors

        Parameters
        ----------
        yls : list
            The y positions of the ray on the left mirror
        yrs : list
            The y positions of the ray on the right mirror
        n : int
            The number of times to propogate the ray

        Returns
        -------
        Ys_left, Ys_right : tuple
            The y positions of the ray on the left and right mirrors
        """
        points_left = np.array([-np.ones(n + 1) * self.L / 2, yls]).T
        points_right = np.array([np.ones(n) * self.L / 2, yrs]).T
        Ys_left = [p[1] for p in points_left]
        Ys_right = [p[1] for p in points_right]
        return Ys_left, Ys_right

    def propogate(
        self,
        pos0,
        n=100,
        show_metadata=True,
        show_fig=True,
        return_ys=False,
        return_fig=False,
        show_start_point=False,
        fig_name=None,
    ):
        """Propogates the ray through the resonator

        Parameters
        ----------
        pos0 : np.array
            The position vector of the ray
        n : int, optional
            The number of times to propogate, by default 1
        show_metadata : bool, optional
            Whether to show the metadata, by default True
        return_ys : bool, optional
            Whether to return the y positions of the ray on the mirrors, by default True
        return_fig : bool, optional
            Whether to return the figure, by default False
        show_start_point : bool, optional
            Whether to show the start point of the ray, by default False
        fig_name : str, optional
            The name of the figure to save, by default None (does not save)

        Returns
        -------
        (Ys_left, Ys_right) or (ax, fig) or None
            The y positions of the ray on the left and right mirrors or the figure and axis or None
        """
        if not self.is_stable:
            print("Resonator is not stable. Plotted ray may not be accurate.")

        self.y0 = pos0[0]
        self.x0 = self.get_left_x(self.y0)
        self.theta0 = pos0[1]

        self.draw(show=False, show_metadata=show_metadata)

        yls, yrs = self.__propogate(pos0, n)
        Ys_left, Ys_right = self.__preprocess(yls, yrs, n)
        Xs_left_calculated = [self.get_left_x(y) for y in Ys_left]
        Xs_right_calculated = [self.get_right_x(y) for y in Ys_right]

        if show_start_point:
            self.ax.plot(self.x0, self.y0, marker="o", ms=10, color="k")
        for i in range(n):
            self.ax.plot(
                (Xs_left_calculated[i], Xs_right_calculated[i]),
                (Ys_left[i], Ys_right[i]),
                c=self.RAY_COLOR,
                lw=0.5,
            )
            self.ax.plot(
                (Xs_right_calculated[i], Xs_left_calculated[i + 1]),
                (Ys_right[i], Ys_left[i + 1]),
                c=self.RAY_COLOR,
                lw=0.5,
            )
        if fig_name is not None:
            plt.savefig(fig_name, dpi=300)
        if show_fig:
            plt.show()
        if return_ys:
            return yls, yrs
        if return_fig:
            return self.fig, self.ax
