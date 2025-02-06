defmodule Wavelets.Analysis do
  @moduledoc """
  Tools for analyzing and visualizing wavelet transforms
  """

  @doc """
  Computes energy distribution in wavelet coefficients or signals
  """
  @spec energy_distribution(list(number) | list(list(number))) :: float | map()
  def energy_distribution(coeffs) when is_list(coeffs) do
    case is_list(hd(coeffs)) do
      true ->
        # 2D case
        coeffs
        |> List.flatten()
        |> compute_energy()
        |> Kernel./(length(List.flatten(coeffs)))

      false ->
        # Compute raw energy
        raw_energy = compute_energy(coeffs)

        # If length <= 2, it's wavelet coefficients
        if length(coeffs) <= 2 do
          # Preserve energy for wavelet coefficients
          raw_energy
        else
          # Energy density for regular signals
          raw_energy / length(coeffs)
        end
    end
  end

  @doc """
  Computes Shannon entropy of wavelet coefficients
  """
  @spec entropy(list(number) | list(list(number))) :: float
  def entropy(coeffs) when is_list(coeffs) do
    coeffs
    |> List.flatten()
    |> Enum.map(&abs/1)
    |> normalize()
    |> Enum.reject(&(&1 <= 0))
    |> Enum.map(fn p -> -p * :math.log2(p) end)
    |> Enum.sum()
  end

  defp compute_energy(signal) do
    signal
    |> Enum.map(&(&1 * &1))
    |> Enum.sum()
  end

  defp normalize(values) do
    sum = values |> Enum.map(&abs/1) |> Enum.sum()

    case sum do
      0.0 -> values
      _ -> Enum.map(values, &(abs(&1) / sum))
    end
  end
end
