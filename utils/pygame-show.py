from numpy._typing import NDArray
import pygame
import numpy as np
import csv
import sys


def open_file(file_name: str) -> tuple[NDArray, int]:
    print(f"Opening file: {file_name}, it might take a while.", end="\r")

    with open(file_name, "r") as f:
        data = list(csv.reader(f, delimiter=","))

    data = np.array([list(map(float, x)) for x in data])

    print("Read completed!                                          ")
    return data, int(np.max(data[:, :1])) + 1


def prepare_coordinates(
    data: NDArray, width: int, height: int
) -> tuple[NDArray, NDArray, NDArray]:
    print("Preparing data.", end="\r")
    # Normalize x, y columns
    data[:, 1:2] *= (width / 2) / np.max(data[:, 1:2])  # * 8/10
    data[:, 2:3] *= (height / 2) / np.max(data[:, 2:3])  # * 8/10

    c_x = width // 2
    c_y = height // 2

    # Cast a intero, conversione a vettore e posizionamento nella schermata
    xs = data[:, 1:2].astype(int).flatten() + c_x
    ys = data[:, 2:3].astype(int).flatten() + c_y

    if np.max(data[:, 3:4]) != 0:
        min_z = np.min(data[:, 3:4])
        if min_z < 0:
            data[:, 3:4] += -min_z
        else:
            data[:, 3:4] -= min_z
        data[:, 3:4] *= 255 / np.max(data[:, 3:4])
        zs = data[:, 3:4].astype(int).flatten()
    else:
        zs = np.full(len(xs), 255)

    print("Data preparation completed!")

    return xs, ys, zs


def play(
    xs: NDArray, ys: NDArray, zs: NDArray, n_corps: int, width: int, height: int
) -> None:
    # Inizializzazione di Pygame
    pygame.init()

    # Creazione della finestra di gioco
    screen = pygame.display.set_mode((width, height))

    # Ciclo principale del gioco
    running = True
    clock = pygame.time.Clock()
    current_frame = 0

    font = pygame.font.SysFont("Arial", 18)

    out_bodies = set()
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        # Pulizia dello schermo
        screen.fill((0, 0, 0))

        frame_text = font.render(
            "Time Step: " + str(current_frame + 0), True, (255, 255, 255)
        )
        screen.blit(frame_text, (10, 10))

        # Aggiornamento delle posizioni dei punti per il frame corrente
        start = n_corps * current_frame
        end = start + n_corps

        if start >= len(xs):
            print("Simulation Completed. Closing.")
            break

        # Disegno dei punti sulla finestra di gioco
        for i in range(start, end):
            x = xs[i]
            y = ys[i]
            if i == start:
                col = (255, 0, 0)
            else:
                # col = (255, 255, 255)
                col = (zs[i], zs[i], zs[i])

            if x > width or y > height:
                n = (i + 1) % n_corps
                if n not in out_bodies:
                    out_bodies.add(n)

                    print( # Se rientra nello schermo e poi riesce non è rilevato
                        f"WARNING: body {n} is off-screen: x = {x}/{width}, y = {y}/{height}",
                        file=sys.stderr,
                    )
                    break
            pygame.draw.circle(screen, col, (x, y), 3, width=0) # circonferenza
            pygame.draw.circle(screen, (0, 0, 0), (x, y), 2, width=0) # interno

        # Aggiornamento del frame corrente
        current_frame += 1

        # Aggiornamento della finestra di gioco
        pygame.display.flip()

        # Limitazione del framerate
        clock.tick(30*10)

    # Uscita dal gioco
    pygame.quit()


def main():
    argc = len(sys.argv)
    if argc < 2 or argc == 3:  # or missing height
        print(
            f"Usage: python3 {sys.argv[0]} FILE_NAME [WIDTH, HEIGHT]", file=sys.stderr
        )
        exit(1)

    if argc == 4:
        width = int(sys.argv[2])
        height = int(sys.argv[3])
    else:
        width = 1024
        height = 768

    data, n_corps = open_file(sys.argv[1])
    xs, ys, zs = prepare_coordinates(data, width, height)
    play(xs, ys, zs, n_corps, width, height)


if __name__ == "__main__":
    main()
