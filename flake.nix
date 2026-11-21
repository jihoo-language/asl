{
  description = "Alexandrov Section Lambda";

  inputs.nixpkgs.url = "github:nixos/nixpkgs/nixos-26.05";

  outputs = { nixpkgs, ... }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];
      forSystems = f: nixpkgs.lib.genAttrs systems (system: f nixpkgs.legacyPackages.${system});
      package = pkgs: pkgs.rustPlatform.buildRustPackage {
        pname = "asl";
        version = "0.1.0";
        src = pkgs.lib.fileset.toSource {
          root = ./.;
          fileset = pkgs.lib.fileset.unions [
            ./Cargo.lock
            ./Cargo.toml
            ./examples
            ./src
            ./tests
          ];
        };
        cargoLock.lockFile = ./Cargo.lock;
        meta = {
          description = "Interpreter for Alexandrov Section Lambda";
          license = pkgs.lib.licenses.asl20;
          mainProgram = "asl";
        };
      };
    in
    {
      packages = forSystems (pkgs: { default = package pkgs; });
      apps = forSystems (pkgs: {
        default = {
          type = "app";
          program = pkgs.lib.getExe (package pkgs);
          meta.description = "Run an ASL program";
        };
      });
      devShells = forSystems (pkgs: {
        default = pkgs.mkShell {
          packages = with pkgs; [ cargo rustc rustfmt clippy ];
        };
      });
      checks = forSystems (pkgs: {
        default = package pkgs;
        relator-brackets = pkgs.runCommand "asl-relator-brackets" {
          nativeBuildInputs = [ (package pkgs) pkgs.gawk ];
        } ''
          set -o pipefail
          awk -f ${./examples/bracket-relator.awk} | asl - > "$out"
          diff -u - "$out" <<'EXPECTED'
          betti(K) = [1, 4, 1]
          {H(K, 2){0} => S.hi, H(K, 2){g0} => S.hi}
          {H(K, 2){0} => S.hi, H(K, 2){g0} => S.hi}
          {H(K, 2){0} => S.lo, H(K, 2){g0} => S.lo}
          EXPECTED
        '';
        flash-bockstein = pkgs.runCommand "asl-flash-bockstein" {
          nativeBuildInputs = [ (package pkgs) pkgs.gawk ];
        } ''
          set -o pipefail
          {
            awk -v m=1 -f ${./examples/flash-bockstein.awk} | asl -
            awk -v m=97 -f ${./examples/flash-bockstein.awk} | asl -
          } > "$out"
          diff -u - "$out" <<'EXPECTED'
          VBock(One, Flash, [0, 1], 1, 1, 2, 1){g0}
          VBock(One, Flash, [0, 1], 97, 1, 2, 97){g0}
          EXPECTED
          for invalid in 0 -1 2.5 bad 1000000; do
            if awk -v m="$invalid" -f ${./examples/flash-bockstein.awk} > rejected.asl; then
              exit 1
            fi
            test ! -s rejected.asl
          done
        '';
        staircase-transfer = pkgs.runCommand "asl-staircase-transfer" {
          nativeBuildInputs = [ (package pkgs) pkgs.gawk ];
        } ''
          set -o pipefail
          {
            awk -v n=3 -f ${./examples/staircase-transfer.awk} | asl -
            awk -v n=6 -f ${./examples/staircase-transfer.awk} | asl -
            awk -v local=1 -v n=7 -f ${./examples/staircase-transfer.awk} | asl -
            awk -v local=1 -v n=20 -f ${./examples/staircase-transfer.awk} | asl -
          } > "$out"
          diff -u - "$out" <<'EXPECTED'
          RVTor(Stair, 4, 6){g0}
          RVTor(Stair, 10, 12){g0}
          RVMTor(Stair, 12, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]){g0}
          RVMTor(Stair, 38, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]){g0}
          EXPECTED
          for invalid in 2 32 2.5 bad; do
            if awk -v n="$invalid" -f ${./examples/staircase-transfer.awk} > rejected.asl; then
              exit 1
            fi
            test ! -s rejected.asl
          done
          awk -v n=20 -v pages=1 -f ${./examples/staircase-transfer.awk} | asl - > bar-pages.out
          test "$(tail -n 1 bar-pages.out)" = 'RVBarE(StairBar, 0, 20, 19, 39, 1){g0}'
          awk -v n=20 -v pages=1 -v dual=1 -v omega=1 -f ${./examples/staircase-transfer.awk} | asl - > bar-cobar.out
          tail -n 3 bar-cobar.out > bar-cobar.last
          {
            echo 'RVBarE(StairBar, 0, 20, 19, 39, 1){g0}'
            echo 'RVPathExt(StairBar, 0, 20, 40){0}'
            tail -n 1 "$out"
          } | diff -u - bar-cobar.last
          if awk -v local=2 -f ${./examples/staircase-transfer.awk} > rejected.asl; then
            exit 1
          fi
          test ! -s rejected.asl
        '';

        sheaf-transgression = pkgs.runCommand "asl-sheaf-transgression" {
          nativeBuildInputs = [ (package pkgs) pkgs.gawk ];
        } ''
          set -o pipefail
          {
            awk -v r=6 -v prime=3 -f ${./examples/sheaf-transgression.awk} | asl -
            awk -v r=6 -v prime=3 -v untwisted=1 -f ${./examples/sheaf-transgression.awk} | asl -
          } > "$out"
          diff -u - "$out" <<'EXPECTED'
          shedim[K, 2, 0, 5] = 1
          ShE(K, 6, 6, 0){2*g0}
          shedim[K, 7, 0, 5] = 0
          shedim[K, 7, 6, 0] = 0
          shedim[K, 2, 0, 5] = 1
          ShE(K, 6, 6, 0){0}
          shedim[K, 7, 0, 5] = 1
          shedim[K, 7, 6, 0] = 1
          EXPECTED
        '';
      });
    };
}
