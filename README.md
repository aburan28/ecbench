# ecbench

The `ecbench` command is copied from [`aburan28/crypto`](https://github.com/aburan28/crypto). Its CLI source lives in `src/bin/`; the benchmark methods, workloads, records, and checkers are imported from the pinned `crypto` library. Changes to those shared APIs belong in `crypto` so there is one implementation.

## Build

```sh
./tools/bootstrap_crypto.sh
cargo build --locked --release
cargo run --locked --release -- methods
```

The bootstrap fetches a shallow, sparse checkout at the SHA in the script into ignored `deps/crypto`. Review that SHA when updating the CLI, and update the copy and pin together. The source reference is [`crypto`'s ecbench documentation](https://github.com/aburan28/crypto/tree/main/docs/ecbench).

The CLI accepts the same specs and session directories as the command in `crypto`. Generated sessions and comparison outputs are evidence; keep their input, receipt, and replay requirements from the source documentation.
