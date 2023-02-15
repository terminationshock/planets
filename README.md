# Planets

With this desktop application, you can simulate the solar system. Observe how the planets move or investigate how a heavy object passing by disturbs the orbits.

![Screenshot](/screenshot.png?raw=true)

## Features

This program
- renders the ecliptic, so that the solar system is seen from above.
- uses accurate physical properties of the planets and their orbits.
- calculates the Keplerian orbits of all planets using Runge-Kutta 4th order integration.
- conserves energy and angular momentum sufficiently well.
- shows orbital parameters as tooltips of the planets.
- handles collisions between planets.

Note that the planet images are not to scale.

## Key bindings

| Key | Description |
| --- | ----------- |
| `+`/`-` | Speed up/slow down the time |
| `Page up`/`Page down` | Zoom in/out |
| Arrow keys | Move |
| `0` | Stop/resume |
| `Home` | Reset view |
| `o` | Toggle orbits |
| `Esc` | Quit |

## Dependencies

- Python: [Pygame](https://www.pygame.org/) and [NumPy](https://numpy.org/)
- [GNU Fortran compiler](https://gcc.gnu.org/fortran/)

## Run

Build the Fortran helper library first with `make`. Then start the application with `./planets.py`.
