# Example Workflows

- `madtq.yaml`: Small workflow with 3 steps

1. `run_madness`: Computes a tequila molecule with madness as backend using the customImage `kottmanj/madness-tequila`
2. `compute-pno-upccd`: Computes a small VQE with `tequila` using the PNO-UpCCD ansatz
3. `convert-to-of`: Computes the (fermionic) Hamiltonian as OpenFermion InteractionOperator and saves as json

-`madoq.yaml`: Small workflow that illustrates the usage in an orquestra workflow similar to [this one](https://github.com/zapatacomputing/z-quantum-vqe/blob/0fdd2a9d7745cd5acab24a7dbe53fb7c0a435e8f/examples/hydrogen-vqe.yaml#L67-L76)
