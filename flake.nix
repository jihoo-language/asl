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
      });
    };
}
