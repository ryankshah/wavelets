defmodule Wavelets.Analysis do
  @moduledoc """
  Tools for analyzing and visualizing wavelet transforms
  """

  @type coefficient_list :: list(number) | list(list(number))
  @type energy_result :: %{level: integer, subband: atom, energy: float} | float

  @doc """
  Computes energy distribution across wavelet coefficients.
  For 1D signals returns a float value representing total energy.
  For 2D signals returns a map with level, subband, and energy information.
  """
  @spec energy_distribution(coefficient_list()) :: energy_result()
  def energy_distribution(coeffs) when is_list(coeffs) do
    cond do
      is_number(hd(coeffs)) ->
        # 1D case
        compute_1d_energy(coeffs)

      is_list(hd(coeffs)) ->
        # 2D or higher case
        compute_2d_energy(coeffs)
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
    |> Enum.map(fn p ->
      if p > 0, do: -p * :math.log2(p), else: 0
    end)
    |> Enum.sum()
  end

  @doc """
  Performs statistical analysis of wavelet coefficients
  """
  @spec statistics(list(number) | list(list(number))) :: map()
  def statistics(coeffs) when is_list(coeffs) do
    flat_coeffs = List.flatten(coeffs)

    %{
      mean: mean(flat_coeffs),
      variance: variance(flat_coeffs),
      skewness: skewness(flat_coeffs),
      kurtosis: kurtosis(flat_coeffs),
      max_amplitude: Enum.max_by(flat_coeffs, &abs/1),
      sparsity: sparsity(flat_coeffs)
    }
  end

  @doc """
  Generates visualization data for plotting
  """
  @spec visualization_data(list(number) | list(list(number))) :: map()
  def visualization_data(coeffs) when is_list(coeffs) do
    cond do
      is_number(hd(coeffs)) ->
        # 1D visualization data
        generate_1d_viz_data(coeffs)

      is_list(hd(coeffs)) ->
        # 2D visualization data
        generate_2d_viz_data(coeffs)
    end
  end

  # Helper functions
  defp compute_1d_energy(coeffs) do
    coeffs
    |> Enum.map(&(&1 * &1))
    |> Enum.sum()
  end

  defp compute_2d_energy(coeffs) do
    coeffs
    |> Enum.map(&compute_1d_energy/1)
    |> Enum.sum()
  end

  defp normalize(values) do
    sum = Enum.sum(Enum.map(values, &abs/1))
    Enum.map(values, &(abs(&1) / sum))
  end

  defp mean(values) do
    Enum.sum(values) / length(values)
  end

  defp variance(values) do
    m = mean(values)

    values
    |> Enum.map(&:math.pow(&1 - m, 2))
    |> mean()
  end

  defp skewness(values) do
    m = mean(values)
    v = variance(values)
    std = :math.sqrt(v)

    values
    |> Enum.map(&:math.pow((&1 - m) / std, 3))
    |> mean()
  end

  defp kurtosis(values) do
    m = mean(values)
    v = variance(values)
    std = :math.sqrt(v)

    values
    |> Enum.map(&:math.pow((&1 - m) / std, 4))
    |> mean()
    # Excess kurtosis
    |> Kernel.-(3)
  end

  defp sparsity(values) do
    threshold = 0.001 * Enum.max_by(values, &abs/1)
    count_small = Enum.count(values, &(abs(&1) < threshold))
    count_small / length(values)
  end

  defp generate_1d_viz_data(coeffs) do
    %{
      values: coeffs,
      x_range: 0..(length(coeffs) - 1) |> Enum.to_list(),
      type: "1d"
    }
  end

  defp generate_2d_viz_data(coeffs) do
    height = length(coeffs)
    width = length(hd(coeffs))

    %{
      values: coeffs,
      dimensions: {width, height},
      type: "2d"
    }
  end
end
