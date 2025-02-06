defmodule Wavelets.MixProject do
  use Mix.Project

  @version "0.1.0"
  @github_url "https://github.com/yourusername/wavelets"

  def project do
    [
      app: :wavelets,
      version: @version,
      elixir: "~> 1.14",
      start_permanent: Mix.env() == :prod,
      deps: deps(),
      description: description(),
      package: package(),
      docs: docs(),
      test_coverage: [tool: ExCoveralls],
      preferred_cli_env: [
        coveralls: :test,
        "coveralls.detail": :test,
        "coveralls.html": :test,
        "coveralls.json": :test,
        "coveralls.post": :test
      ],
      dialyzer: [
        plt_add_apps: [:mix, :ex_unit],
        plt_file: {:no_warn, "priv/plts/dialyzer.plt"}
      ],
      elixirc_paths: elixirc_paths(Mix.env())
    ]
  end

  def application do
    [
      extra_applications: [:logger]
    ]
  end

  defp deps do
    [
      # Development dependencies
      {:ex_doc, "~> 0.30", only: :dev, runtime: false},
      {:dialyxir, "~> 1.4", only: [:dev, :test], runtime: false},
      {:credo, "~> 1.7", only: [:dev, :test], runtime: false},
      {:nx, "~> 0.6.2"},

      # Test dependencies
      {:excoveralls, "~> 0.18", only: :test},
      {:stream_data, "~> 0.6", only: [:dev, :test]}
    ]
  end

  defp description do
    """
    A comprehensive wavelet transform library for Elixir, providing implementations
    of various wavelet transforms including DWT, CWT, and wavelet packets, with
    support for both built-in and custom wavelet filters.
    """
  end

  defp package do
    [
      name: "wavelets",
      files: ~w(lib .formatter.exs mix.exs README* LICENSE*),
      licenses: ["MIT"],
      links: %{"GitHub" => @github_url}
    ]
  end

  defp docs do
    [
      main: "readme",
      source_url: @github_url,
      source_ref: "v#{@version}",
      extras: ["README.md", "LICENSE.md"],
      groups_for_modules: [
        "Core Transforms": [
          Wavelets.DWT,
          Wavelets.IDWT,
          Wavelets.NDWT,
          Wavelets.SWT,
          Wavelets.CWT,
          Wavelets.WaveletPacket
        ],
        "Wavelet Filters": [
          Wavelets.Filter,
          Wavelets.Filters.Daubechies,
          Wavelets.Filters.Coiflets,
          Wavelets.Filters.Symlet,
          Wavelets.Filters.Biorthogonal,
          Wavelets.Filters.Morlet
        ],
        "Analysis Tools": [
          Wavelets.Analysis,
          Wavelets.WaveletFunctions
        ],
        Utilities: [
          Wavelets.Utils.Complex,
          Wavelets.Utils.Precision
        ]
      ]
    ]
  end

  defp elixirc_paths(:test), do: ["lib", "test/support"]
  defp elixirc_paths(_), do: ["lib"]
end
