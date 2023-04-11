# CaLiG

## Run

After compiling cpp file `calig.cpp` to `calig`

```shell
g++ -O2 ./calig.cpp -o ./calig -std=c++11
```

run

```shell
./calig -d <data_path> -q <query_path> -s <stream_path>
```

- <data_path>: use the path of data graph to replace it.
- <query_path>: use the path of query graph to replace it.
- <stream_path>: use the path of stream graph to replace it.

## Example

```shell
./calig -d ./data/test_G -q ./data/test_Q -s ./data/test_S
```

## Format

### data graph / query graph

```
v <id> <label>
...
e <id1> <id2>
...
```

- Lines starting with "v" represent vertices;
- Lines starting with "e" represent edges.

### stream

```
e <id1> <id2>
...
e <-id1-1> <-id2-1>
...
```

- `e <id1> <id2>` represents the addition of the edge `(id1, id2)`;
- `e <-id1-1> <-id2-1>` represents the deletion of the edge `(id1, id2)`.